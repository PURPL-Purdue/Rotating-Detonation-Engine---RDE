function paths = MoC_Trace_Characteristics(field, geometry)
% MoC_Trace_Characteristics Extracts and clips characteristic line segments (C- and C+)
% from the solved Method of Characteristics field for visualization.
%
% Bounded by the chamber top exit (y <= yMax) and the downstream shock line
% instead of being abruptly cut at xMax.
%
% Inputs:
%   field: Struct containing solved field points (.X, .Y, .Parent_C_Minus, .Parent_C_Plus)
%   geometry: Struct containing geometric boundaries (.xMax, .chamber_length or .yMax, .y_inj, .beta/.shock_angle)
%
% Output:
%   paths: Struct containing:
%          .X, .Y: NaN-delimited coordinates for all characteristic segments within the domain
%          .X_minus, .Y_minus: C- family segments
%          .X_plus, .Y_plus:   C+ family segments

    if isfield(geometry, 'xMax')
        xMax = geometry.xMax;
    else
        xMax = 140.44;
    end

    if isfield(geometry, 'chamber_length')
        yMax = geometry.chamber_length;
    elseif isfield(geometry, 'yMax')
        yMax = geometry.yMax;
    else
        yMax = 29.21;
    end

    if isfield(geometry, 'y_inj')
        y_inj = geometry.y_inj;
    elseif isfield(geometry, 'y_slip0')
        y_inj = geometry.y_slip0;
    else
        y_inj = 7.37;
    end

    if isfield(geometry, 'shock_angle') && geometry.shock_angle > 50
        beta = geometry.shock_angle;
    elseif isfield(geometry, 'beta')
        beta = geometry.beta;
    else
        beta = 70.33;
    end

    num_pts = length(field.X);
    
    % Pre-allocate coordinate accumulators (3 elements per segment: [p1, p2, NaN])
    X_minus = nan(1, num_pts * 3);
    Y_minus = nan(1, num_pts * 3);
    minus_idx = 0;

    X_plus = nan(1, num_pts * 3);
    Y_plus = nan(1, num_pts * 3);
    plus_idx = 0;

    for i = 1:num_pts
        % C- Characteristic: from Parent_C_Minus to current point i
        pm = field.Parent_C_Minus(i);
        if ~isnan(pm) && pm > 0 && pm <= num_pts
            x1 = field.X(pm); y1 = field.Y(pm);
            x2 = field.X(i);  y2 = field.Y(i);
            [xc, yc, vis] = MoC_Clip_Line(x1, y1, x2, y2, xMax, yMax, y_inj, beta);
            if vis
                X_minus(minus_idx + (1:3)) = [xc(1), xc(2), NaN];
                Y_minus(minus_idx + (1:3)) = [yc(1), yc(2), NaN];
                minus_idx = minus_idx + 3;
            end
        end

        % C+ Characteristic: from Parent_C_Plus to current point i
        pp = field.Parent_C_Plus(i);
        if ~isnan(pp) && pp > 0 && pp <= num_pts
            x1 = field.X(pp); y1 = field.Y(pp);
            x2 = field.X(i);  y2 = field.Y(i);
            [xc, yc, vis] = MoC_Clip_Line(x1, y1, x2, y2, xMax, yMax, y_inj, beta);
            if vis
                X_plus(plus_idx + (1:3)) = [xc(1), xc(2), NaN];
                Y_plus(plus_idx + (1:3)) = [yc(1), yc(2), NaN];
                plus_idx = plus_idx + 3;
            end
        end
    end

    % Trim to actual segment count
    paths.X_minus = X_minus(1:minus_idx);
    paths.Y_minus = Y_minus(1:minus_idx);
    paths.X_plus  = X_plus(1:plus_idx);
    paths.Y_plus  = Y_plus(1:plus_idx);

    % Combined paths
    paths.X = [paths.X_minus, paths.X_plus];
    paths.Y = [paths.Y_minus, paths.Y_plus];
end
