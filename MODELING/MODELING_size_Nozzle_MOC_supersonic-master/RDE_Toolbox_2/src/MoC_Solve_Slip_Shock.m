function ShockSlipField = MoC_Solve_Slip_Shock(ShockSeeds, geom, field, gas, options)
% Steady irrotational MoC in the same wave-fixed frame as MoC_Solve_Field.
% Shock and slip remain straight. Each interior point has C+ and C- parents.
% One incoming C- row starts at each accepted shock seed and ends at the slip.
%
% K+ = Theta - Nu and K- = Theta + Nu are constant on their Mach lines.
% Angles are degrees throughout, matching MoC_Solve_Field.
% After constructing the net, P0 and T0 are mapped from the shock along
% assumed straight streamlines PARALLEL TO THE SLIP. Static properties vary
% with the local MoC Mach number. These mapping lines are not characteristics.
%
% This is the reduced-model approximation described by Grunenwald (2023),
%
% Inputs: gas.gamma, gas.R, gas.P0, gas.T0; totals refer to the upstream gas.
% Output: ordinary point arrays, with local parent indices and ShockSeed IDs.
% X,Y are unwrapped input coordinates (mm). Theta is degrees; P is Pa; T is K.
%
% Approximation: a prescribed shock need not satisfy outgoing compatibility,
% and a prescribed slip need not match product pressure. Diagnostics report
% those mismatches. Skipped shock seeds leave an approximate gap in the input
% data; their original indices and reasons are returned in Diagnostics.
% Below the first supplied streamline, use the first accepted seed's totals.
%
% options.maxSteps limits the number of shock-seeded characteristic rows.
% There is no fixed stepSize: intersections determine the point spacing.

    %% 1. Solver Settings and Point Storage 
    % This initializes all the limits, geometry, output arrays and diagnostics. 
    if nargin < 5
        options = struct();
    end
    if ~isfield(options, 'maxSteps')
        options.maxSteps = 100000;
    end
    if isfield(options,'stepSize') && ~isempty(options.stepSize)
        warning('MoC:SlipShockStepSizeIgnored', ...
            'stepSize is unused: characteristic intersections now set the spacing.');
    end

    max_steps = options.maxSteps; % limits the number of shock seeds
    MAX_POINTS = 1000000;         % limits the number of points in the shock slip region

    validateattributes(max_steps, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(gas.gamma, {'numeric'}, {'scalar', 'finite', '>', 1});
    validateattributes(gas.R, {'numeric'}, {'scalar', 'finite', 'positive'});
    validateattributes(gas.P0, {'numeric'}, {'scalar', 'finite', 'positive'});
    validateattributes(gas.T0, {'numeric'}, {'scalar', 'finite', 'positive'});

    gamma = gas.gamma;
    slip_angle = geom.slip_angle;
    shock_angle = geom.shock_angle; % Already includes phi.
    assert(0 < slip_angle && slip_angle < shock_angle && shock_angle < 90, ...
        'MoC:SlipShockInput', 'Require 0 < slip angle < shock angle < 90 degrees.');
    assert(geom.yMax > geom.y_inj, 'MoC:SlipShockInput', 'Exit must be above the triple point.');

    % COORDINATE SHIFTING: 

    % Rotate coordinates so that they align with the slip line,
    %  so the straight slip becomes march_y = 0.
    % This is not a change of reference frame or an additional flow rotation.
    % This transformation is applied later around lines 90-100
    % This makes the slip boundary essentially act as an y axis
    % creates new x', y', and theta' coordinates
    slip_cos = cosd(slip_angle);
    slip_sin = sind(slip_angle);
    shock_slope = tand(shock_angle - slip_angle);
    end_x = (geom.yMax - geom.y_inj) / slip_sin; 
    position_tolerance = 1e-10 * max(1, end_x);

    point_names = {'X','Y','M','Theta','Nu','P','T','rho','P0','T0', ...
        'a','V','v_x','v_y','Region','Type','Parent_C_Minus','Parent_C_Plus','ShockSeed'};
    capacity = max(10000, numel(ShockSeeds.X));
    for k = 1:numel(point_names)
        ShockSlipField.(point_names{k}) = nan(capacity, 1);
    end
    idx = 1;

    diagnostics.completed = false;
    diagnostics.status = 'No shock seeds.';
    diagnostics.steps = 0;
    diagnostics.startStation = NaN;
    diagnostics.lastStation = NaN;
    diagnostics.maxShockCompatibilityDeg = 0;
    diagnostics.maxSlipPressureRelativeMismatch = NaN;
    diagnostics.slipComparisonCount = 0;
    diagnostics.slipBoundaryPoints = 0;
    diagnostics.reflectedPlusSamples = 0;
    diagnostics.maxSlipCompatibilityDeg = 0;
    diagnostics.skippedShockSeeds = [];
    diagnostics.skippedSeedReasons = {};
    diagnostics.failedRows = 0;
    diagnostics.totalConditionExtrapolationCount = 0;

    %% 2. Load the Supplied Shock Seeds

    seed_count = numel(ShockSeeds.X);

    assert(seed_count <= MAX_POINTS, 'MoC:SlipShockInput', 'Too many shock seeds.');
    for name = {'X','Y','M','Theta','P0_Ratio'}

        values = ShockSeeds.(name{1});

            assert(numel(values) == seed_count && isreal(values) && all(isfinite(values(:))), ...
                'MoC:SlipShockInput', 'Shock seed arrays must be finite, real, and equal length.');
    end

    % Transforming the shock seeds to match the local coordinate system created above.
    seed_x = (ShockSeeds.X(:) - geom.xMax)*slip_cos + ...
             (ShockSeeds.Y(:) - geom.y_inj)*slip_sin;

    seed_y = -(ShockSeeds.X(:) - geom.xMax)*slip_sin + ...
              (ShockSeeds.Y(:) - geom.y_inj)*slip_cos;


            assert(all(seed_x >= -position_tolerance) && ...
                all(abs(seed_y - shock_slope*seed_x) < 1e-6*max(1,end_x)), ...
                'MoC:SlipShockInput', 'Seeds must lie on the unwrapped straight shock.');

    
    % Adding the supplied shock seeds from the MoC_Solve_Field function

    % sorting the seeds by x'
    [seed_x, order] = sort(seed_x);
    seed_y = seed_y(order);
    seed_mach = ShockSeeds.M(order);
    seed_mach = seed_mach(:);

    % Transforming flow angle where theta' = theta - slip angle
    seed_state.Theta = ShockSeeds.Theta(order);
    seed_state.Theta = seed_state.Theta(:) - slip_angle;
    seed_state.M = seed_mach;
    seed_state.Nu = pmAngle(max(seed_mach,1),gamma);

    % Importing and computing shock thermodynamics for the seed points.
    % Rankine-Hugionot relation has already been computed in MoC_Oblique_Shock_Jump.m
    seed_state.P0 = gas.P0 * ShockSeeds.P0_Ratio(order);
    seed_state.P0 = seed_state.P0(:);
    seed_state.T0 = repmat(gas.T0, seed_count, 1);

    if isfield(ShockSeeds, 'P0')
        seed_state.P0 = ShockSeeds.P0(order);
        seed_state.P0 = seed_state.P0(:);
    end
    if isfield(ShockSeeds, 'T0')
        seed_state.T0 = ShockSeeds.T0(order);
        seed_state.T0 = seed_state.T0(:);
    end
    assert(all(isfinite(seed_state.P0) & seed_state.P0 > 0) && ...
        all(isfinite(seed_state.T0) & seed_state.T0 > 0), ...
        'MoC:SlipShockInput', 'Seed total conditions must be positive and finite.');

    % Merge repeated geometry only when its flow states agree.
    keep = true(seed_count, 1);
    for k = 2:seed_count
        if seed_x(k) - seed_x(k-1) <= position_tolerance
            assert(abs(seed_state.Theta(k)-seed_state.Theta(k-1)) < 1e-8 && ...
                abs(seed_state.M(k)-seed_state.M(k-1)) < 1e-8 && ...
                abs(log(seed_state.P0(k)/seed_state.P0(k-1))) < 1e-8 && ...
                abs(log(seed_state.T0(k)/seed_state.T0(k-1))) < 1e-8, ...
                'MoC:SlipShockInput', 'Conflicting states at duplicate shock seeds.');
            keep(k) = false;
        end
    end
    seed_x = seed_x(keep);
    seed_y = seed_y(keep);
    seed_mach = seed_mach(keep);
    order = order(keep);
    seed_state = selectState(seed_state, find(keep));
    seed_count = numel(seed_x);


    %% 3. Build the Characteristic Net, One Incoming C- Row at a Time
    % Some seeds have two incoming families; that is allowed here. Only
    % discard states that cannot emit a downstream C- into this region.
    mu = asind(1./max(seed_mach,1));
    c_minus_angle = seed_state.Theta - mu;
    usable = seed_mach > 1 & c_minus_angle > -90 + 1e-8 & c_minus_angle < 90 - 1e-8 & ...
        tand(c_minus_angle) < shock_slope;
    diagnostics.skippedShockSeeds = order(~usable);
    reasons = repmat({'C- points out of the shock-slip region.'},seed_count,1);
    reasons(c_minus_angle <= -90 + 1e-8 | c_minus_angle >= 90 - 1e-8) = {'Characteristic runs upstream.'};
    reasons(seed_mach <= 1) = {'Subsonic or sonic shock state.'};
    diagnostics.skippedSeedReasons = reasons(~usable);

    seed_x = seed_x(usable);
    seed_y = seed_y(usable);
    order = order(usable);
    seed_state = selectState(seed_state,find(usable));
    seed_count = numel(seed_x);
    current_front = [];
    stopped = false;

    if seed_count > 0
        % Store the original seeds once. ShockSeed retains the INPUT index,
        % so main can merge accepted roots without duplicating shock points.
        for k = 1:seed_count
            point = selectState(seed_state,k);
            point.x = seed_x(k);
            point.y = seed_y(k);
            storePoint(point,5,NaN,NaN,order(k));
        end
        diagnostics.startStation = seed_x(1);
        row_count = nnz(seed_x <= end_x+position_tolerance);

        for step = 1:min(row_count,max_steps)
            % Like Solve_Field, a front contains point IDs. Here its order
            % follows one C- from the shock down to the slip.
            next_front = zeros(1,numel(current_front)+2);
            next_front(1) = step; % This row's supplied shock seed.
            front_count = 1;
            row_failed = false;
            reached_end = false;

            %% A. Interior Intersections
            for k = 1:numel(current_front)
                p_plus = current_front(k);
                p_minus = next_front(front_count);
                plus = readPoint(p_plus);
                minus = readPoint(p_minus);
                [point,solved] = interiorPoint(plus,minus,gamma);

                % A C+ that already left through the prescribed shock cannot
                % meet this later incoming C- inside the domain.
                outside = point.x <= max(plus.x,minus.x)+position_tolerance || ...
                    point.y < -position_tolerance || ...
                    point.y > shock_slope*point.x+position_tolerance;
                if outside
                    if step > 1
                        residual = shockMismatch(plus,readPoint(step-1), ...
                            readPoint(step),shock_slope);
                        diagnostics.maxShockCompatibilityDeg = max( ...
                            diagnostics.maxShockCompatibilityDeg,residual);
                    end
                    continue
                end
                if ~solved
                    row_failed = true;
                    break
                end
                if point.x > end_x
                    reached_end = true;
                    break
                end
                if idx > MAX_POINTS
                    stopped = true;
                    break
                end

                front_count = front_count+1;
                next_front(front_count) = storePoint(point,0,p_minus,p_plus,order(step));
                if ShockSlipField.Type(p_plus) == 2
                    diagnostics.reflectedPlusSamples = diagnostics.reflectedPlusSamples+1;
                end
            end

            %% B. Slip Boundary and Reflected C+
            % The last C- meets y=0. Tangency supplies theta=0 and incoming
            % compatibility supplies Nu and Mach. Its stored boundary point
            % emits a C+ into the next row; it is not a new arbitrary seed.
            p_minus = next_front(front_count);
            incident = readPoint(p_minus);
            if incident.y > position_tolerance && ~row_failed && ~reached_end && ~stopped
                [point,solved,residual] = slipPoint(incident,gamma);
                if ~solved
                    row_failed = true;
                elseif point.x <= end_x && idx <= MAX_POINTS
                    front_count = front_count+1;
                    next_front(front_count) = storePoint(point,2,p_minus,NaN,order(step));
                    diagnostics.slipBoundaryPoints = diagnostics.slipBoundaryPoints+1;
                    diagnostics.maxSlipCompatibilityDeg = max( ...
                        diagnostics.maxSlipCompatibilityDeg,residual);
                elseif idx > MAX_POINTS
                    stopped = true;
                end
            end

            %% C. Advance to the Next Supplied Shock Seed
            current_front = next_front(1:front_count);
            diagnostics.steps = step;
            diagnostics.lastStation = seed_x(step);
            diagnostics.failedRows = diagnostics.failedRows+row_failed;
            if stopped
                break
            end
        end

        diagnostics.completed = diagnostics.steps == row_count && row_count > 0 && ...
            seed_x(end) >= end_x-position_tolerance && ...
            diagnostics.failedRows == 0 && ~stopped;
        if stopped
            diagnostics.status = 'Maximum point count reached.';
        elseif diagnostics.steps < row_count
            diagnostics.status = 'Maximum characteristic rows reached (maxSteps).';
        elseif diagnostics.failedRows > 0
            diagnostics.status = 'Some characteristic rows could not be continued; returning the solved net.';
        elseif seed_x(end) < end_x-position_tolerance || row_count == 0
            diagnostics.status = 'Shock seed coverage ended before the slip/exit intersection.';
        else
            diagnostics.status = 'Finished the characteristic intersection net.';
        end
    else
        diagnostics.status = 'No usable supersonic shock seeds.';
    end
    if ~isempty(diagnostics.skippedShockSeeds)
        diagnostics.status = sprintf('%s Skipped %d shock seeds; gaps are approximate.', ...
            diagnostics.status,numel(diagnostics.skippedShockSeeds));
        warning('MoC:SlipShockSkippedSeeds','%s',diagnostics.status);
    end

    %% 4. Map Total Conditions and Calculate Static Properties
    for k = 1:numel(point_names)
        ShockSlipField.(point_names{k}) = ShockSlipField.(point_names{k})(1:idx-1);
    end

    % Straight lines parallel to the slip are local_y = constant. Their
    % intersections with the shock local_y = shock_slope*local_x give the
    % source station directly. Map from the ORIGINAL accepted seed totals,
    % not from adjacent characteristic points (which lie on other streamlines).
    local_y = -(ShockSlipField.X-geom.xMax)*slip_sin + ...
               (ShockSlipField.Y-geom.y_inj)*slip_cos;
    source_x = local_y/shock_slope;
    if seed_count > 0
        diagnostics.totalConditionExtrapolationCount = nnz( ...
            source_x < seed_x(1)-position_tolerance | ...
            source_x > seed_x(end)+position_tolerance);
        source_x = min(max(source_x,seed_x(1)),seed_x(end));
        if seed_count == 1
            ShockSlipField.P0(:) = seed_state.P0(1);
            ShockSlipField.T0(:) = seed_state.T0(1);
        else
            ShockSlipField.P0 = interp1(seed_x,seed_state.P0,source_x,'linear');
            ShockSlipField.T0 = interp1(seed_x,seed_state.T0,source_x,'linear');
        end
        % Preserve the supplied shock states exactly.
        ShockSlipField.P0(1:seed_count) = seed_state.P0;
        ShockSlipField.T0(1:seed_count) = seed_state.T0;
    end
    ShockSlipField.P = ShockSlipField.P0 ./ ...
        (1+(gamma-1)*ShockSlipField.M.^2/2).^(gamma/(gamma-1));
    ShockSlipField.T = ShockSlipField.T0 ./ (1+(gamma-1)*ShockSlipField.M.^2/2);
    ShockSlipField.rho = ShockSlipField.P ./ (gas.R*ShockSlipField.T);
    ShockSlipField.a = sqrt(gamma*gas.R*ShockSlipField.T);
    ShockSlipField.V = ShockSlipField.M.*ShockSlipField.a;
    ShockSlipField.v_x = ShockSlipField.V.*cosd(ShockSlipField.Theta);
    ShockSlipField.v_y = ShockSlipField.V.*sind(ShockSlipField.Theta);
    ShockSlipField.InChamber = ShockSlipField.Y <= geom.yMax+position_tolerance;

    %% 5. Compare Pressure Across the Prescribed Slip
    if isfield(field,'Type') && isfield(field,'P') && isfield(geom,'x_slip0')

        product_slip = find(field.Type(:)==2 & isfinite(field.P(:)) & field.P(:)>0);

        [product_x,unique_ids] = unique(field.X(product_slip));

        if numel(product_x) >= 2
            product_pressure = field.P(product_slip(unique_ids));
            slip_ids = find(ShockSlipField.Type==2 & ShockSlipField.InChamber);
            sample_x = ShockSlipField.X(slip_ids)-geom.xMax+geom.x_slip0;
            matched_pressure = interp1(product_x,product_pressure,sample_x,'linear',NaN);
            valid = isfinite(matched_pressure);


            diagnostics.slipComparisonCount = nnz(valid);
            if any(valid)
                diagnostics.maxSlipPressureRelativeMismatch = max(abs( ...
                    ShockSlipField.P(slip_ids(valid))-matched_pressure(valid))./matched_pressure(valid));
            end
        end
    end
    ShockSlipField.Diagnostics = diagnostics;
    if ~diagnostics.completed
        warning('MoC:SlipShockIncomplete','%s',diagnostics.status);
    end
    if diagnostics.maxShockCompatibilityDeg>0.5 || diagnostics.maxSlipPressureRelativeMismatch>0.05
        warning('MoC:SlipShockBoundaryMismatch', ...
            'Fixed boundaries: shock residual %.3g deg; slip pressure mismatch %.3g.', ...
            diagnostics.maxShockCompatibilityDeg,diagnostics.maxSlipPressureRelativeMismatch);
    end


    function point = readPoint(id)
        % Local coordinates are parallel/normal to the straight slip.
        point.x = (ShockSlipField.X(id)-geom.xMax)*slip_cos + ...
                  (ShockSlipField.Y(id)-geom.y_inj)*slip_sin;

        point.y = -(ShockSlipField.X(id)-geom.xMax)*slip_sin + ...
                   (ShockSlipField.Y(id)-geom.y_inj)*slip_cos;

        point.Theta = ShockSlipField.Theta(id)-slip_angle;

        point.M = ShockSlipField.M(id);
        
        point.Nu = ShockSlipField.Nu(id);
    end

    function id = storePoint(point,point_type,parent_minus,parent_plus,seed_id)
        % Grow storage in blocks, then trim once at the end.
        if idx > capacity
            new_capacity = min(MAX_POINTS,2*capacity);
            for n = 1:numel(point_names)
                ShockSlipField.(point_names{n})(capacity+1:new_capacity,1) = NaN;
            end
            capacity = new_capacity;
        end
        id = idx;
        ShockSlipField.X(id) = geom.xMax+point.x*slip_cos-point.y*slip_sin;
        ShockSlipField.Y(id) = geom.y_inj+point.x*slip_sin+point.y*slip_cos;
        ShockSlipField.Theta(id) = point.Theta+slip_angle;
        ShockSlipField.M(id) = point.M;
        ShockSlipField.Nu(id) = point.Nu;
        ShockSlipField.Type(id) = point_type;
        ShockSlipField.Region(id) = 3; % Slip-shock; product points use Region 2.
        ShockSlipField.Parent_C_Minus(id) = parent_minus;
        ShockSlipField.Parent_C_Plus(id) = parent_plus;
        ShockSlipField.ShockSeed(id) = seed_id;
        idx = idx+1;
    end
end

function [point,solved] = interiorPoint(plus,minus,gamma)
% Irrotational compatibility, followed by the intersection of the two rays.
    K_plus = plus.Theta-plus.Nu;
    K_minus = minus.Theta+minus.Nu;
    point.Theta = (K_minus+K_plus)/2;
    point.Nu = (K_minus-K_plus)/2;
    point.M = machFromNu(point.Nu,gamma);
    point.x = NaN;
    point.y = NaN;
    solved = isfinite(point.M) && abs(point.Theta)+asind(1/point.M) < 90-1e-8;
    if ~solved
        return
    end

    % Average endpoint slopes improve geometry without iterating the state.
    mu = asind(1/point.M);
    m_plus = (tand(plus.Theta+asind(1/plus.M))+tand(point.Theta+mu))/2;
    m_minus = (tand(minus.Theta-asind(1/minus.M))+tand(point.Theta-mu))/2;
    point.x = (minus.y-plus.y+m_plus*plus.x-m_minus*minus.x)/(m_plus-m_minus);
    point.y = plus.y+m_plus*(point.x-plus.x);
    solved = isfinite(point.x) && isfinite(point.y);
end

function [point,solved,residual] = slipPoint(incident,gamma)
% At the straight slip, theta=0 and K- = theta+nu comes from the incident ray.
    point.Theta = 0;
    point.Nu = incident.Theta+incident.Nu;
    point.M = machFromNu(point.Nu,gamma);
    point.x = NaN;
    point.y = 0;
    residual = NaN;
    solved = isfinite(point.M);
    if ~solved
        return
    end

    m_minus = (tand(incident.Theta-asind(1/incident.M))-tand(asind(1/point.M)))/2;
    point.x = incident.x-incident.y/m_minus;
    solved = isfinite(point.x) && point.x > incident.x;
    residual = abs(point.Theta+point.Nu-incident.Theta-incident.Nu);
end

function residual = shockMismatch(plus,first,last,shock_slope)
% Compare outgoing K+ with the imposed shock data where the ray exits.
    m_plus = tand(plus.Theta+asind(1/plus.M));
    x = (m_plus*plus.x-plus.y)/(m_plus-shock_slope);
    residual = 0;
    if ~isfinite(x) || x <= plus.x || x < first.x || x > last.x
        return
    end
    fraction = (x-first.x)/(last.x-first.x);
    boundary_K = (1-fraction)*(first.Theta-first.Nu) + ...
        fraction*(last.Theta-last.Nu);
    residual = abs(boundary_K-(plus.Theta-plus.Nu));
end

function M = machFromNu(nu,gamma)
% Reuse the main solver's inverse, with a bracketed fallback near its limits.
% A sonic or impossible PM angle ends this row instead of being clamped.
    nu_limit = 90*(sqrt((gamma+1)/(gamma-1))-1);
    M = NaN;
    if ~isfinite(nu) || nu <= 0 || nu >= nu_limit
        return
    end
    M = MoC_PM_Inv(nu,gamma);
    if ~isreal(M) || ~isfinite(M) || M <= 1 || abs(pmAngle(M,gamma)-nu) > 1e-8
        upper = 2;
        while pmAngle(upper,gamma) < nu
            upper = 2*upper;
        end
        M = fzero(@(mach) pmAngle(mach,gamma)-nu,[1,upper]);
    end
end

function nu = pmAngle(M,gamma)
% Prandtl-Meyer angle in degrees.
    factor = sqrt((gamma+1)/(gamma-1));
    root = sqrt(M.^2-1);
    nu = factor*atand(root/factor)-atand(root);
end

function selected = selectState(state,indices)
    selected.Theta = state.Theta(indices);
    selected.M = state.M(indices);
    selected.Nu = state.Nu(indices);
    selected.P0 = state.P0(indices);
    selected.T0 = state.T0(indices);
end
