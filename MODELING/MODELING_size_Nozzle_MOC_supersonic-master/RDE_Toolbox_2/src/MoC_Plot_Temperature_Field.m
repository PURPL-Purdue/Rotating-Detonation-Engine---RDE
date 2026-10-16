function tempFigure = MoC_Plot_Temperature_Field(field, geometry, xMax, yMax)
% MOC_PLOT_TEMPERATURE_FIELD Plots the solved characteristic mesh points
% colored by local temperature to produce a point-gradient flow visualization.
%
% Inputs:
%   field: Structure containing solved coordinates and properties (.X, .Y, .T)
%   geometry: Structure containing domain parameters (xMax, yMax, y_inj, x_inj, phi, etc.)
%   xMax: (Optional) Chamber circumference [mm]
%   yMax: (Optional) Chamber length [mm]
%
% Output:
%   tempFigure: Handle to the generated figure

    %% 1. Geometric Parameters
    if nargin < 3 || isempty(xMax)
        xMax = geometry.xMax;
    end
    if nargin < 4 || isempty(yMax)
        yMax = geometry.chamber_length;
    end
    
    y_inj = geometry.y_inj;
    x_inj = geometry.x_inj;
    phi = geometry.phi;
    beta = geometry.shock_angle;
    
    if isfield(geometry, 'slip_angle')
        slip_angle = geometry.slip_angle;
    else
        slip_angle = beta;
    end

    y_plot_top = yMax + 10;
    x_infill = xMax - x_inj;
    x_det_floor = y_inj * tand(phi);

    %% 2. Setup Figure and Axes
    % Stagger window position slightly so it doesn't completely occlude the net figure
    tempFigure = figure('Name', 'RDE Temperature Field (Point Gradient)', ...
        'NumberTitle', 'off', ...
        'Position', [120, 160, 1200, 360], ...
        'Color', 'w');
    
    ax = axes('Parent', tempFigure);
    hold(ax, 'on');

    %% 3. Plot Solid Temperature Dots
    if ~isempty(field) && isfield(field, 'T') && isfield(field, 'X') && isfield(field, 'Y')
        % Filter out NaNs
        valid = ~isnan(field.X) & ~isnan(field.Y) & ~isnan(field.T);
        X = field.X(valid);
        Y = field.Y(valid);
        T = field.T(valid);
        
        %% Filter Points to Physical Combustion Chamber Domain
        % 1. Must be below chamber exit (the black line: Y <= yMax)
        % 2. Must not be inside the unburned reactant infill triangle:
        %    The infill hypotenuse is x = x_infill + y / tand(phi) for y <= y_inj
        in_infill_c1 = (Y <= y_inj + 1e-5) & (X >= (x_infill + Y ./ tand(phi) - 1e-4));
        x_shock_bnd_c1 = xMax + max(0, Y - y_inj) / tand(beta);
        in_c1 = (Y <= yMax + 1e-4) & ~in_infill_c1 & (X <= x_shock_bnd_c1 + 1e-4);
        
        X_c1 = X(in_c1);
        Y_c1 = Y(in_c1);
        T_c1 = T(in_c1);

        % Tunable dot size (increase for denser appearance, decrease if overlapping excessively)
        dot_size = 14; 
        
        % Plot Cycle 1 (Primary domain: 0 to xMax, clipped to physical boundaries)
        scatter(ax, X_c1, Y_c1, dot_size, T_c1, 'filled', ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 1.0); % Solid opaque dots
        
        % Plot Cycle 2 (Periodic domain: shifted by xMax, shown up to 1.5 * xMax)
        X2 = X + xMax;
        Y2 = Y;
        in_infill_c2 = (Y2 <= y_inj + 1e-5) & (X2 >= (x_infill + xMax + Y2 ./ tand(phi) - 1e-4));
        in_c2 = (Y2 <= yMax + 1e-4) & (X2 <= 1.5 * xMax + 1e-3) & ~in_infill_c2;
        
        if any(in_c2)
            scatter(ax, X2(in_c2), Y2(in_c2), dot_size, T(in_c2), 'filled', ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', 1.0);
        end
        
        % Colorbar / Temperature Key
        cb = colorbar(ax);
        ylabel(cb, 'Static Temperature [K]', 'FontSize', 10, 'FontWeight', 'bold');
        
        % Colormap: 'turbo' gives a rich, continuous thermal gradient
        colormap(ax, 'turbo');
    end

    %% 4. Overlay Domain Boundaries for Context
    % Infill Triangle Shaded Region (Unburned Fresh Reactants)
    patch(ax, [x_infill, xMax, xMax], [0, y_inj, 0], [0.95, 0.75, 0.75], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.6);
    
    % Detonation Wave (Cycle 1)
    plot(ax, [x_det_floor, 0], [0, y_inj], 'k-', 'LineWidth', 1.5);
    
    % Infill Triangle Boundary (Cycle 1)
    plot(ax, [x_infill, xMax], [0, y_inj], 'm-', 'LineWidth', 2.0);
    plot(ax, [xMax, xMax], [0, y_inj], 'k-', 'LineWidth', 1.0);

    % Oblique Shock and Slip Lines (Cycle 1)
    x_shock_0 = (y_plot_top - y_inj) / tand(beta);
    x_slip_0  = (y_plot_top - y_inj) / tand(slip_angle);
    plot(ax, [0, x_shock_0], [y_inj, y_plot_top], 'b-', 'LineWidth', 1.2);
    plot(ax, [0, x_slip_0],  [y_inj, y_plot_top], 'b--', 'LineWidth', 1.2);

    % Cycle 2 Detonation and Shock Boundaries (shifted by xMax)
    plot(ax, [x_det_floor + xMax, xMax], [0, y_inj], 'k-', 'LineWidth', 1.5);
    plot(ax, [xMax, xMax + x_shock_0], [y_inj, y_plot_top], 'b-', 'LineWidth', 1.2);
    plot(ax, [xMax, xMax + x_slip_0],  [y_inj, y_plot_top], 'b--', 'LineWidth', 1.2);

    % Chamber Exit Boundary
    plot(ax, [0, 1.5 * xMax], [yMax, yMax], 'k-', 'LineWidth', 2.0);

    %% 5. Formatting and Aspect Ratio
    xlim(ax, [0, 1.5 * xMax]);
    ylim(ax, [0, y_plot_top]);
    yticks(ax, 0:10:y_plot_top);
    pbaspect(ax, [1.5 * xMax, y_plot_top, 1]);
    
    grid(ax, 'on');
    box(ax, 'on');
    
    xlabel(ax, 'Axial Distance (Unwrapped) [mm]');
    ylabel(ax, 'Height [mm]');
    title(ax, 'RDE Temperature Field [K] (Point Gradient)', 'FontWeight', 'bold');
end

