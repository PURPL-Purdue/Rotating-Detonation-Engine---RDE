function [fieldFigure, netFigure] = MoC_Plot_Field(field, geometry, xMax, yMax)
% MoC_Plot_Field Plots the unwrapped Method of Characteristics net and RDE flow domain.
%
% Inputs:
%   field: Structure containing the solved characteristic points and properties
%   geometry: Structure containing domain parameters (xMax, yMax, y_inj, x_inj, phi, etc.)
%   xMax: (Optional) Unwrapped chamber circumference [mm]
%   yMax: (Optional) Chamber length [mm]
%
% Outputs:
%   fieldFigure: Handle to the generated figure
%   netFigure: Handle to the generated figure

    %% 1. Extract Geometric Parameters
    xMax = geometry.xMax;
    yMax = geometry.chamber_length;
    y_inj = geometry.y_inj;
    x_inj = geometry.x_inj;
    phi = geometry.phi;
    beta = geometry.shock_angle;
    slip_angle = geometry.slip_angle;

    %% 2. Setup Figure and Axes
    fig = figure('Name', 'Method of Characteristics Net', 'NumberTitle', 'off', ...
        'Position', [100, 200, 1200, 320]);
    ax = axes('Parent', fig);
    hold(ax, 'on');

    x_infill = xMax - x_inj;
    x_det_floor = y_inj * tand(phi);
    y_plot_top = yMax + 10;

    %% 3. Infill Triangle (Cycle 1)
    % Shaded coral/light-red region
    patch(ax, [x_infill, xMax, xMax], [0, y_inj, 0], [1, 0.5, 0.5], 'EdgeColor', 'none');
    
    % Infill boundary (Magenta)
    plot(ax, [x_infill, xMax], [0, y_inj], 'm', 'LineWidth', 2.5);
    
    % Infill rear vertical boundary / wave boundary (Black)
    plot(ax, [xMax, xMax], [0, y_inj], 'k', 'LineWidth', 1.5);

    %% 4. Detonation Waves (Cycle 1 and Cycle 2)
    % Cycle 1 detonation line at x = 0
    plot(ax, [x_det_floor, 0], [0, y_inj], 'k', 'LineWidth', 1);

    %% 5. Oblique Shock and Slip Lines (Blue)
    % Cycle 1: from triple point (0, y_inj)
    x_shock_0 = (y_plot_top - y_inj) / tand(beta);
    x_slip_0  = (y_plot_top - y_inj) / tand(slip_angle);
    plot(ax, [0, x_shock_0], [y_inj, y_plot_top], 'Color', [0 0 1], 'LineWidth', 1.8);
    plot(ax, [0, x_slip_0],  [y_inj, y_plot_top], 'Color', [0 0 1], 'LineWidth', 1.8);

    % Cycle 2: from triple point (xMax, y_inj)
    plot(ax, [xMax, xMax + x_shock_0], [y_inj, y_plot_top], 'Color', [0 0 1], 'LineWidth', 1.8);
    plot(ax, [xMax, xMax + x_slip_0],  [y_inj, y_plot_top], 'Color', [0 0 1], 'LineWidth', 1.8);

    %% 6. Chamber Exit Boundary (Thick Horizontal Black Line)
    plot(ax, [0, 1.5 * xMax], [yMax, yMax], 'k', 'LineWidth', 2);

    %% 7. Characteristic Net Lines
    if ~isempty(field) && isfield(field, 'X') && ~isempty(field.X)
        % Ensure geometry passed to tracer contains chamber boundaries
        geom_trace = geometry;
        geom_trace.xMax = xMax;
        geom_trace.chamber_length = yMax;
        geom_trace.y_inj = y_inj;
        geom_trace.shock_angle = beta;
        
        paths = MoC_Trace_Characteristics(field, geom_trace);
        
        % Gray characteristic net lines
        netColor = [0.8 0.8 0.8];
        
        % Plot Cycle 1 characteristic net
        plot(ax, paths.X, paths.Y, 'Color', netColor, 'LineWidth', 0.5);
        
        % Plot Cycle 2 characteristic net (shifted by 1 full circumference xMax)
        plot(ax, paths.X + xMax, paths.Y, 'Color', netColor, 'LineWidth', 0.5);
    end

    %% 8. Axes Formatting and 1:1 Aspect Ratio Scaling
    xlim(ax, [0, 1.5 * xMax]);
    ylim(ax, [0, y_plot_top]);
    yticks(ax, 0:10:30);
    pbaspect(ax, [1.5 * xMax, y_plot_top, 1]);
    
    grid(ax, 'on');
    box(ax, 'on');
    
    xlabel(ax, 'Axial Distance (Unwrapped)');
    ylabel(ax, 'Height');
    title(ax, 'Method of Characteristics Net', 'FontWeight', 'bold');

    fieldFigure = fig;
    netFigure = fig;
end
