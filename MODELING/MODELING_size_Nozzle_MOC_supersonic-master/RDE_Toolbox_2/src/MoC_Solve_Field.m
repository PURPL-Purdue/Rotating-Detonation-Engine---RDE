function [Field, ShockSeeds] = MoC_Solve_Field(seedPoints, geom, gas)
% MOC_Solve_Field Solves the 2D Method of Characteristics supersonic flow field
% Uses a "Working Front" 1D-Array architecture to handle complex boundaries,
% free-jet infill regions, and open engine tops without complex i,j indexing.
% 
% Essentially runs a C- characteristic line from the bottom-most point in the 
% working front, and traces it until it hits a boundary, then runs a C+ from 
% the same point and a C- from the next point above it, and determines that 
% interaction, then rinse and repeat for all the points and finally, 
% run a single c+ from the topmost point.

    %% 1. Pre-allocate 1D arrays for better performance and set Piecewise Functions
    MAX_PTS = 10000000;
    
    Field.X = nan(MAX_PTS, 1);
    Field.Y = nan(MAX_PTS, 1);
    Field.M = nan(MAX_PTS, 1);
    Field.Theta = nan(MAX_PTS, 1); % Flow angle [degrees]
    Field.Nu = nan(MAX_PTS, 1); % Prandtl-Meyer angle [degrees]
    Field.P = nan(MAX_PTS, 1); % Static pressure [Pa]
    Field.T = nan(MAX_PTS, 1); % Static temperature [K]
    Field.rho = nan(MAX_PTS, 1); % Density [kg/m^3]
    Field.Region = nan(MAX_PTS, 1); % Region for plotting purposes: 1=Post-Det, 2=Slip-Shock, 3=Out of Top
    
    % Point types: 0=Interior, 1=Floor, 2=Slipline, 3=Infill, 4=Engine Top, 5=Shock
    Field.Type = zeros(MAX_PTS, 1); 
    
    % Genealogy for plotting the characteristic net later
    Field.Parent_C_Minus = nan(MAX_PTS, 1);
    Field.Parent_C_Plus  = nan(MAX_PTS, 1);

    % Physical chamber exit height (boundary between Region 1 and Region 3)
    if isfield(geom, 'chamber_length')
        y_exit = geom.chamber_length;
    elseif isfield(geom, 'yMax')
        y_exit = geom.yMax;
    else
        y_exit = 29.21;
    end
    
    % Triple point height
    if isfield(geom, 'y_inj')
        y_inj = geom.y_inj;
    elseif isfield(geom, 'y_slip0')
        y_inj = geom.y_slip0;
    else
        y_inj = 5.5;
    end
    
    % Shock wave angle and slope
    if isfield(geom, 'shock_angle')
        beta_shock = geom.shock_angle;
    else
        beta_shock = 70.33;
    end
    
    if isfield(geom, 'phi')
        phi = geom.phi;
    elseif isfield(geom, 'infill_phi')
        phi = geom.infill_phi;
    else
        phi = 0;
    end
    
    m_shock = tand(beta_shock + phi);

    % Piecewise functions for the shock and slip lines (valid for x >= geom.xMax)
    shock_line = @(x) (y_inj + m_shock * (x - geom.xMax)) .* (x >= geom.xMax);
    slip_line  = @(x) (y_inj + tand(geom.slip_angle + phi) * (x - geom.xMax)) .* (x >= geom.xMax);

    % Pre-allocate Shock Seeds struct for Region 2 evaluation
    MAX_SHOCK_SEEDS = 10000;
    shock_idx = 0;
    ShockSeeds.X = nan(MAX_SHOCK_SEEDS, 1);
    ShockSeeds.Y = nan(MAX_SHOCK_SEEDS, 1);
    ShockSeeds.M = nan(MAX_SHOCK_SEEDS, 1);
    ShockSeeds.Theta = nan(MAX_SHOCK_SEEDS, 1);
    ShockSeeds.Nu = nan(MAX_SHOCK_SEEDS, 1);
    ShockSeeds.P0_Ratio = nan(MAX_SHOCK_SEEDS, 1);
    ShockSeeds.Parent_Ray = nan(MAX_SHOCK_SEEDS, 1);

    %% 2. Load the Initial Value Line (IVL)
    % Assuming seedPoints contains arrays X, Y, M, Theta, Nu ordered from bottom to top
    num_seeds = length(seedPoints.X);
    idx = 1;
    
    for i = 1:num_seeds
        Field.X(idx)     = seedPoints.X(i);
        Field.Y(idx)     = seedPoints.Y(i);
        Field.M(idx)     = seedPoints.M(i);
        Field.Theta(idx) = seedPoints.Theta(i);
        Field.Nu(idx)    = seedPoints.Nu(i);
        
        % Identify if the seed point is already on a bottom boundary
        if seedPoints.Y(i) <= 1e-6
            if seedPoints.X(i) < geom.x_infill_start
                Field.Type(idx) = 1; % Floor
            else
                Field.Type(idx) = 3; % Infill
            end
        else
            Field.Type(idx) = 0; % Interior
        end
        Field.Region(idx) = 1; % Seed points are in Post-Det zone (Region 1)
        
        idx = idx + 1;
    end
    
    % The working front is an array of indices spanning from bottom to top of the chamber
    current_front = 1:num_seeds; 
    last_infill_idx = NaN; % Tracks the last point generated on the infill boundary

    %% 3. Marching Loop (Loop until the whole front passes the next triple point)
    while min(Field.X(current_front)) < geom.xMax
        
        % If all lines exited the top of the engine, the front is empty.
        if isempty(current_front)
            fprintf('FRONT IS EMPTY at min X = %.2f\n', min(Field.X(1:idx-1)));
            break; 
        end
        
        next_front = [];
        
        %% A. Bottom Boundary (Floor or Infill)
        % If the lowest point on the current front is an INTERIOR point, it must 
        % emit a C- ray downwards to hit a bottom boundary.
        p_bottom = current_front(1); % Bottom point of the current working front
        
        if Field.Type(p_bottom) == 0 % Case if the bottom point is an internal point (this should almost always happen)
            
            % Extract properties of the downward C- ray
            theta1 = Field.Theta(p_bottom);
            nu1    = Field.Nu(p_bottom);
            mu1    = asind(1 / Field.M(p_bottom)); % Mach angle of the downward C- ray (relative to theta)
            m_minus_slope = tand(theta1 - mu1); % Slope of the downward C- characteristic line (relative to horizontal)
            
            x1 = Field.X(p_bottom);
            y1 = Field.Y(p_bottom);
            
            % Where does this ray hit y=0?
            x_zero_crossing = x1 - y1 / m_minus_slope;
            
            if x_zero_crossing < geom.x_infill_start % If the C- ray hits the solid engine floor before the initial infill boundary
                % --- HIT THE SOLID ENGINE FLOOR ---
                x_int = x_zero_crossing;
                y_int = 0;
                
                if x_int >= x1 - 1e-5 % Ensure the x-intersection is not behind the starting point
                    theta_int = 0; % Solid wall boundary condition
                    
                    % Invariant along C-: J_minus = Theta + Nu
                    J_minus = theta1 + nu1;
                    nu_int  = J_minus - theta_int;
                    
                    % Add to Field
                    Field.X(idx) = x_int; Field.Y(idx) = y_int;
                    Field.Theta(idx) = theta_int; Field.Nu(idx) = nu_int;
                    Field.M(idx) = MoC_PM_Inv(nu_int, gas.gamma); % Helper func needed
                    Field.Type(idx) = 1; % Floor
                    Field.Region(idx) = 1;
                    Field.Parent_C_Minus(idx) = p_bottom;
                    
                    next_front = [next_front, idx];
                    idx = idx + 1;
                end
                
            else
                % --- HIT THE RIGID WEDGE INFILL BOUNDARY ---
                % Boundary condition: Rigid wedge -> Constant Theta
                theta_int = geom.infill_phi; % Angle of the injection triangle
                
                % Invariant along C-: J_minus = Theta + Nu
                J_minus = theta1 + nu1;
                nu_int = J_minus - theta_int;
                
                m_jet = tand(theta_int);
                
                % Get the previous infill point to draw the boundary line
                if isnan(last_infill_idx)
                    x_prev = geom.x_infill_start;
                    y_prev = 0;
                else
                    x_prev = Field.X(last_infill_idx);
                    y_prev = Field.Y(last_infill_idx);
                end
                
                % Intersect C- line with Jet Boundary line
                x_int = (m_minus_slope*x1 - y1 - m_jet*x_prev + y_prev) / (m_minus_slope - m_jet);
                y_int = m_jet * (x_int - x_prev) + y_prev;
                
                if x_int >= x1 - 1e-5
                    % Add to Field Struct
                    Field.X(idx) = x_int; 
                    Field.Y(idx) = y_int;
                    Field.Theta(idx) = theta_int; 
                    Field.Nu(idx) = nu_int;
                    Field.M(idx) = MoC_PM_Inv(nu_int, gas.gamma);
                    Field.Type(idx) = 3; % Infill
                    Field.Region(idx) = 1;
                    Field.Parent_C_Minus(idx) = p_bottom;
                    
                    last_infill_idx = idx; % Save for the next infill ray
                    next_front = [next_front, idx];
                    idx = idx + 1;
                end
            end
        end
        
        %% B. Interior Intersections
        % Loop through adjacent pairs on the current front
        for k = 1:(length(current_front)-1)
            p_up   = current_front(k);   % Emits C+ (upward), k represents the current point in the front, 
                                         % k=1 only has a C+ ray because the C- associated with it has 
                                         % already been drawn above
            p_down = current_front(k+1); % Emits C- (downward)
            
            % C+ Properties
            theta_up = Field.Theta(p_up); 
            nu_up = Field.Nu(p_up);
            mu_up = asind(1 / Field.M(p_up));
            m_plus = tand(theta_up + mu_up);
            
            % C- Properties
            theta_dn = Field.Theta(p_down); 
            nu_dn = Field.Nu(p_down);
            mu_dn = asind(1 / Field.M(p_down));
            m_minus_slope = tand(theta_dn - mu_dn);
            
            % Intersect lines
            x1 = Field.X(p_up); 
            y1 = Field.Y(p_up);
            x2 = Field.X(p_down); 
            y2 = Field.Y(p_down);
            
            x_int = (m_plus*x1 - y1 - m_minus_slope*x2 + y2) / (m_plus - m_minus_slope);
            y_int = m_plus * (x_int - x1) + y1;
            
            % SHOCK HANDLING: If resultant point is behind either of its parent points, then delete it because it is invalid
            if x_int < x1 - 1e-5 || x_int < x2 - 1e-5
                continue; % Skip adding this point to next_front!
            end
            
            % Calculate interior flow properties using invariants
            J_plus  = theta_up - nu_up;
            J_minus = theta_dn + nu_dn;
            
            theta_int = (J_minus + J_plus) / 2;
            nu_int    = (J_minus - J_plus) / 2;
            
            if nu_int < 0.1
                nu_int = 0.1; % Cap to prevent subsonic breakdown, and prevent the return of imaginary Mach numbers and errors
            end
            
            % Only keep the point if it hasn't escaped the top of the analysis zone
            if y_int <= geom.y_engine_top
                % --- CHECK IF INTERIOR INTERSECTION CROSSED THE SHOCK LINE ---
                % Shock exists for x >= geom.xMax and y >= y_inj
                hit_shock = false;
                if y_int >= y_inj - 1e-5 && x_int >= geom.xMax - 1e-5
                    x_shock_at_y = geom.xMax + (y_int - y_inj) / m_shock;
                    if x_int >= x_shock_at_y - 1e-4
                        hit_shock = true;
                    end
                end
                
                if hit_shock
                    % Upward C+ ray from p_up intersects the shock line
                    denom = m_plus - m_shock;
                    if abs(denom) > 1e-6
                        x_s = (m_plus * x1 - y1 - m_shock * geom.xMax + y_inj) / denom;
                        y_s = y_inj + m_shock * (x_s - geom.xMax);
                        
                        if x_s >= geom.xMax - 1e-4 && y_s >= y_inj - 1e-4 && y_s <= geom.y_engine_top
                            % Rankine-Hugoniot jump across oblique shock
                            M_up = Field.M(p_up);
                            theta_up_val = Field.Theta(p_up);
                            [M2, theta2, nu2, shockProps] = MoC_Oblique_Shock_Jump(M_up, theta_up_val, beta_shock + phi, gas.gamma);
                            
                            % Store shock boundary point in Field
                            Field.X(idx) = x_s;
                            Field.Y(idx) = y_s;
                            Field.Theta(idx) = theta2;
                            Field.Nu(idx) = nu2;
                            Field.M(idx) = M2;
                            Field.Type(idx) = 5; % Shock
                            Field.Region(idx) = 2; % Region 2: Slip-Shock boundary
                            Field.Parent_C_Plus(idx) = p_up;
                            Field.Parent_C_Minus(idx) = NaN;
                            
                            % Store in ShockSeeds accumulator
                            shock_idx = shock_idx + 1;
                            if shock_idx <= MAX_SHOCK_SEEDS
                                ShockSeeds.X(shock_idx) = x_s;
                                ShockSeeds.Y(shock_idx) = y_s;
                                ShockSeeds.M(shock_idx) = M2;
                                ShockSeeds.Theta(shock_idx) = theta2;
                                ShockSeeds.Nu(shock_idx) = nu2;
                                ShockSeeds.P0_Ratio(shock_idx) = shockProps.P02_P01;
                                ShockSeeds.Parent_Ray(shock_idx) = p_up;
                            end
                            idx = idx + 1;
                        end
                    end
                    % Ray terminated at shock boundary; do not add to next_front
                    continue;
                end

                Field.X(idx) = x_int; 
                Field.Y(idx) = y_int;
                Field.Theta(idx) = theta_int; 
                Field.Nu(idx) = nu_int;
                Field.M(idx) = MoC_PM_Inv(nu_int, gas.gamma);
                Field.Type(idx) = 0; % Interior
                
                % Region classification: Region 1 (Chamber) vs Region 3 (Above exit)
                if y_int <= geom.yMax
                    Field.Region(idx) = 1; % Post-Detonation chamber zone
                else
                    Field.Region(idx) = 3; % Out of top of engine
                end
                
                Field.Parent_C_Plus(idx) = p_up;
                Field.Parent_C_Minus(idx) = p_down;
                
                next_front = [next_front, idx];
                idx = idx + 1;
            else
                % Ray from p_up goes out the top of the analysis zone
                x_top = x1 + (geom.y_engine_top - y1) / m_plus;
                Field.X(idx) = x_top; Field.Y(idx) = geom.y_engine_top;
                Field.Theta(idx) = theta_int; Field.Nu(idx) = nu_int;
                Field.M(idx) = MoC_PM_Inv(nu_int, gas.gamma);
                Field.Type(idx) = 4; % Engine Top 
                Field.Region(idx) = 3; % Out of top
                Field.Parent_C_Plus(idx) = p_up;
                Field.Parent_C_Minus(idx) = NaN; % No C- parent
                % Do not add to next_front as it exits the analysis zone
                idx = idx + 1;
            end
        end
        
        %% C. Top Boundary (Slipline)
        % If the highest point is an INTERIOR point, it must emit a C+ ray upwards.
        p_top = current_front(end);
        
        if Field.Type(p_top) == 0 % Interior
            theta1 = Field.Theta(p_top);
            nu1    = Field.Nu(p_top);
            mu1    = asind(1 / Field.M(p_top));
            m_plus = tand(theta1 + mu1);
            
            x1 = Field.X(p_top); 
            y1 = Field.Y(p_top);
            
            % Slipline is an inviscid solid wall at a specific angle
            theta_slip = geom.slip_angle;
            m_slip = tand(theta_slip);
            
            % Slipline geometrical definition: y - y_slip0 = m_slip * (x - x_slip0)
            x_int = (m_plus*x1 - y1 - m_slip*geom.x_slip0 + geom.y_slip0) / (m_plus - m_slip);
            y_int = m_plus * (x_int - x1) + y1;
            
            % Check if top ray crosses the shock before/at the slip line
            hit_shock_top = false;
            if y_int >= y_inj - 1e-5 && x_int >= geom.xMax - 1e-5
                x_shock_at_y = geom.xMax + (y_int - y_inj) / m_shock;
                if x_int >= x_shock_at_y - 1e-4
                    hit_shock_top = true;
                end
            end
            
            if hit_shock_top
                denom = m_plus - m_shock;
                if abs(denom) > 1e-6
                    x_s = (m_plus * x1 - y1 - m_shock * geom.xMax + y_inj) / denom;
                    y_s = y_inj + m_shock * (x_s - geom.xMax);
                    
                    if x_s >= geom.xMax - 1e-4 && y_s >= y_inj - 1e-4 && y_s <= geom.y_engine_top
                        M_up = Field.M(p_top);
                        theta_up_val = Field.Theta(p_top);
                        [M2, theta2, nu2, shockProps] = MoC_Oblique_Shock_Jump(M_up, theta_up_val, beta_shock + phi, gas.gamma);
                        
                        Field.X(idx) = x_s;
                        Field.Y(idx) = y_s;
                        Field.Theta(idx) = theta2;
                        Field.Nu(idx) = nu2;
                        Field.M(idx) = M2;
                        Field.Type(idx) = 5; % Shock
                        Field.Region(idx) = 2; % Region 2: Slip-Shock boundary
                        Field.Parent_C_Plus(idx) = p_top;
                        Field.Parent_C_Minus(idx) = NaN;
                        
                        shock_idx = shock_idx + 1;
                        if shock_idx <= MAX_SHOCK_SEEDS
                            ShockSeeds.X(shock_idx) = x_s;
                            ShockSeeds.Y(shock_idx) = y_s;
                            ShockSeeds.M(shock_idx) = M2;
                            ShockSeeds.Theta(shock_idx) = theta2;
                            ShockSeeds.Nu(shock_idx) = nu2;
                            ShockSeeds.P0_Ratio(shock_idx) = shockProps.P02_P01;
                            ShockSeeds.Parent_Ray(shock_idx) = p_top;
                        end
                        idx = idx + 1;
                    end
                end
                % Top ray terminated at shock; do not add to next_front
            elseif y_int <= geom.y_engine_top
                % Calculate properties
                theta_int = theta_slip;
                J_plus = theta1 - nu1;
                nu_int = theta_int - J_plus;
                
                Field.X(idx) = x_int; 
                Field.Y(idx) = y_int;
                Field.Theta(idx) = theta_int; 
                Field.Nu(idx) = nu_int;
                Field.M(idx) = MoC_PM_Inv(nu_int, gas.gamma);
                Field.Type(idx) = 2; % Slipline
                
                if y_int <= geom.yMax
                    Field.Region(idx) = 1; % Post-Detonation chamber zone
                else
                    Field.Region(idx) = 3; % Out of top
                end
                
                Field.Parent_C_Plus(idx) = p_top;
                
                next_front = [next_front, idx];
                idx = idx + 1;
            else
                % Ray from p_top goes out the top
                x_top = x1 + (geom.y_engine_top - y1) / m_plus;
                Field.X(idx) = x_top; Field.Y(idx) = geom.y_engine_top;
                Field.Theta(idx) = theta1; 
                Field.Nu(idx) = nu1; % Just carry properties
                Field.M(idx) = MoC_PM_Inv(nu1, gas.gamma);
                Field.Type(idx) = 4; % Engine Top
                Field.Region(idx) = 3; % Out of top
                Field.Parent_C_Plus(idx) = p_top;
                idx = idx + 1;
            end
        end
        
        %% D. Advance the front
        current_front = next_front;
        
    end
    
    %% 4. Trim the pre-allocated arrays
    Field.X = Field.X(1:idx-1);
    Field.Y = Field.Y(1:idx-1);
    Field.M = Field.M(1:idx-1);
    Field.Theta = Field.Theta(1:idx-1);
    Field.Nu = Field.Nu(1:idx-1);
    Field.Type = Field.Type(1:idx-1);
    Field.Region = Field.Region(1:idx-1);
    Field.Parent_C_Minus = Field.Parent_C_Minus(1:idx-1);
    Field.Parent_C_Plus = Field.Parent_C_Plus(1:idx-1);

    % Trim ShockSeeds and store inside Field
    if shock_idx > 0
        ShockSeeds.X = ShockSeeds.X(1:shock_idx);
        ShockSeeds.Y = ShockSeeds.Y(1:shock_idx);
        ShockSeeds.M = ShockSeeds.M(1:shock_idx);
        ShockSeeds.Theta = ShockSeeds.Theta(1:shock_idx);
        ShockSeeds.Nu = ShockSeeds.Nu(1:shock_idx);
        ShockSeeds.P0_Ratio = ShockSeeds.P0_Ratio(1:shock_idx);
        ShockSeeds.Parent_Ray = ShockSeeds.Parent_Ray(1:shock_idx);
    else
        ShockSeeds.X = [];
        ShockSeeds.Y = [];
        ShockSeeds.M = [];
        ShockSeeds.Theta = [];
        ShockSeeds.Nu = [];
        ShockSeeds.P0_Ratio = [];
        ShockSeeds.Parent_Ray = [];
    end
    Field.ShockSeeds = ShockSeeds;

end
