function [] = MoC_Plot_Field(field, geometry, plots_requested, plot_options)
% MoC_Plot_Field Plots the unwrapped Method of Characteristics net and RDE flow domain.
%
% Inputs:
%   field: Structure containing the solved characteristic points and properties
%   geometry: Structure containing domain parameters (xMax, yMax, y_inj, x_inj, phi, etc.)
%   plots_requested: Cell array of strings specifying which plots to generate (e.g., {'field', 'net'})
%       1 = Method of Characteristics Net
%       2 = Static Temperature Field [K]
%       3 = Static Pressure Field [bar]
%       4 = Mach Number Field
%       5 = Flow Angle (Theta) Field [deg]
%       6 = Density Field [kg/m^3]
%       7 = Velocity Magnitude Field [m/s]

%   plot_options: Vector of numbers specifying which plot options to include
%       1 = Show Slip Line
%       2 = Show Shock Line
%       3 = Show Borders
%       4 = Show Infill Region
%       5 = Show Detonation Triangle

%% Extract parameters

Y = field.Y;
X = field.X;
Region = field.Region;
plot_logic = field.plot_logic; % Logical index for points within the chamber (Region 1 and 2)
datapoint_indices = field.datapoint_indices; % Index numbers for all points in the field
y_inj = geometry.y_inj;
phi = geometry.phi;
beta = geometry.beta;
delta = geometry.delta;
yMax = geometry.yMax;
xMax = geometry.xMax;
chamber_length = geometry.chamber_length;
x_infill_start = geometry.x_infill_start;

%% Plot 

for i = 1:length(plots_requested)
    switch plots_requested(i)
        case 1

            base_figure = figure('Name', 'Method of Characteristics Net', 'NumberTitle', 'off', ...
            'Position', [100, 200, 1200, 320]);
            ax = axes('Parent', base_figure);
            hold(ax, 'on');

            idx = find(plot_logic);
            h_pts = plot(ax, X(idx), Y(idx), 'k.', 'MarkerSize', 5);
            if ~isempty(Region)
                idx_reg2 = idx(Region(idx) == 2);
                plot(ax, X(idx_reg2), Y(idx_reg2), 'b.', 'MarkerSize', 10);
            end

            % Label points with their index number (for coarse/medium meshes <= 3000 pts)
            if length(idx) <= 3000
                labels = " " + string(datapoint_indices(idx));
                text(ax, X(idx), Y(idx), labels, ...
                    'FontSize', 7, ...
                    'Color', [0.25 0.25 0.25], ...
                    'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment', 'left');
            end

            % Enable data tip index on click/hover for any mesh size
            row_idx = dataTipTextRow('Point Index', datapoint_indices(idx));
            h_pts.DataTipTemplate.DataTipRows(end+1) = row_idx;

            plot([0, (yMax - y_inj)/tand(beta + phi)], [y_inj, yMax], 'r--', 'LineWidth', 2); % Shock line
            plot([0, (yMax - y_inj)/tand(delta + phi)], [y_inj, yMax], 'g--', 'LineWidth', 2); % Slip line
            patch([x_infill_start, xMax, xMax, x_infill_start], [0, 0, y_inj, 0], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            patch([0, y_inj*tand(phi), 0, 0], [y_inj, 0, 0, y_inj], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle

            xlim(ax, [0, xMax]);
            ylim(ax, [0, chamber_length]);
            yticks(ax, 0:10:30);
            pbaspect(ax, [xMax, chamber_length + 10, 1]);
                
            grid(ax, 'on');
            box(ax, 'on');
                
            xlabel(ax, 'Axial Distance - mm (Unwrapped)');
            ylabel(ax, 'Height - mm');
            title(ax, ['Method of Characteristics Net - ' num2str(length(X)) ' Points'], 'FontWeight', 'bold');
            
        case 2
            temps = field.T;

            Static_Temperature_Field = figure('Name', 'Static Temperature Field', 'NumberTitle', 'off', ...
            'Position', [100, 200, 1200, 320]);
            ax = axes('Parent', Static_Temperature_Field);
            hold(ax, 'on');
            
            scatter(ax, X(plot_logic), Y(plot_logic), 20, temps(plot_logic), 'filled');
            c = colorbar;               % Displays the color scale bar
            c.Label.String = 'Temperature (K)'; % Labels the color key axis
            colormap('turbo'); 

            plot([0, (yMax - y_inj)/tand(beta + phi)], [y_inj, yMax], 'r--', 'LineWidth', 2); % Shock line
            plot([0, (yMax - y_inj)/tand(delta + phi)], [y_inj, yMax], 'g--', 'LineWidth', 2); % Slip line
            patch([x_infill_start, xMax, xMax, x_infill_start], [0, 0, y_inj, 0], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            patch([0, y_inj*tand(phi), 0, 0], [y_inj, 0, 0, y_inj], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle

            xlim(ax, [0, xMax]);
            ylim(ax, [0, chamber_length]);
            yticks(ax, 0:10:30);
            pbaspect(ax, [xMax, chamber_length + 10, 1]);
                
            grid(ax, 'on');
            box(ax, 'on');
                
            xlabel(ax, 'Axial Distance - mm (Unwrapped)');
            ylabel(ax, 'Height - mm');
            title(ax, ['Static Temperature Field - ' num2str(sum(plot_logic)) ' Points'], 'FontWeight', 'bold');

        case 3
            pressures = field.P;

            Static_Pressure_Field = figure('Name', 'Static Pressure Field', 'NumberTitle', 'off', ...
            'Position', [100, 200, 1200, 320]);
            ax = axes('Parent', Static_Pressure_Field);
            hold(ax, 'on');
            
            scatter(ax, X(plot_logic), Y(plot_logic), 20, pressures(plot_logic) / 1e5, 'filled');
            c = colorbar;               % Displays the color scale bar
            c.Label.String = 'Pressure (bar)'; % Labels the color key axis
            colormap('turbo'); 

            plot([0, (yMax - y_inj)/tand(beta + phi)], [y_inj, yMax], 'r--', 'LineWidth', 2); % Shock line
            plot([0, (yMax - y_inj)/tand(delta + phi)], [y_inj, yMax], 'g--', 'LineWidth', 2); % Slip line
            patch([x_infill_start, xMax, xMax, x_infill_start], [0, 0, y_inj, 0], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            patch([0, y_inj*tand(phi), 0, 0], [y_inj, 0, 0, y_inj], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle

            xlim(ax, [0, xMax]);
            ylim(ax, [0, chamber_length]);
            yticks(ax, 0:10:30);
            pbaspect(ax, [xMax, chamber_length + 10, 1]);
                
            grid(ax, 'on');
            box(ax, 'on');
                
            xlabel(ax, 'Axial Distance - mm (Unwrapped)');
            ylabel(ax, 'Height - mm');
            title(ax, ['Static Pressure Field - ' num2str(sum(plot_logic)) ' Points'], 'FontWeight', 'bold');

        case 4

            mach = field.M;

            Static_Mach_Field = figure('Name', 'Mach Number Field', 'NumberTitle', 'off', ...
            'Position', [100, 200, 1200, 320]);
            ax = axes('Parent', Static_Mach_Field);
            hold(ax, 'on');
            
            scatter(ax, X(plot_logic), Y(plot_logic), 20, mach(plot_logic), 'filled');
            c = colorbar;               % Displays the color scale bar
            c.Label.String = 'Mach Number'; % Labels the color key axis
            colormap('turbo'); 

            plot([0, (yMax - y_inj)/tand(beta + phi)], [y_inj, yMax], 'r--', 'LineWidth', 2); % Shock line
            plot([0, (yMax - y_inj)/tand(delta + phi)], [y_inj, yMax], 'g--', 'LineWidth', 2); % Slip line
            patch([x_infill_start, xMax, xMax, x_infill_start], [0, 0, y_inj, 0], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            patch([0, y_inj*tand(phi), 0, 0], [y_inj, 0, 0, y_inj], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle

            xlim(ax, [0, xMax]);
            ylim(ax, [0, chamber_length]);
            yticks(ax, 0:10:30);
            pbaspect(ax, [xMax, chamber_length + 10, 1]);
                
            grid(ax, 'on');
            box(ax, 'on');
                
            xlabel(ax, 'Axial Distance - mm (Unwrapped)');
            ylabel(ax, 'Height - mm');
            title(ax, ['Mach Number Field - ' num2str(sum(plot_logic)) ' Points'], 'FontWeight', 'bold');

        case 5

        flow_angle = field.Theta;

            Flow_Angle = figure('Name', 'Flow Angle Field', 'NumberTitle', 'off', ...
            'Position', [100, 200, 1200, 320]);
            ax = axes('Parent', Flow_Angle);
            hold(ax, 'on');
            
            scatter(ax, X(plot_logic), Y(plot_logic), 20, flow_angle(plot_logic), 'filled');
            c = colorbar;               % Displays the color scale bar
            c.Label.String = 'Flow Angle'; % Labels the color key axis
            colormap('turbo'); 

            plot([0, (yMax - y_inj)/tand(beta + phi)], [y_inj, yMax], 'r--', 'LineWidth', 2); % Shock line
            plot([0, (yMax - y_inj)/tand(delta + phi)], [y_inj, yMax], 'g--', 'LineWidth', 2); % Slip line
            patch([x_infill_start, xMax, xMax, x_infill_start], [0, 0, y_inj, 0], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            patch([0, y_inj*tand(phi), 0, 0], [y_inj, 0, 0, y_inj], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle

            xlim(ax, [0, xMax]);
            ylim(ax, [0, chamber_length]);
            yticks(ax, 0:10:30);
            pbaspect(ax, [xMax, chamber_length + 10, 1]);
                
            grid(ax, 'on');
            box(ax, 'on');
                
            xlabel(ax, 'Axial Distance - mm (Unwrapped)');
            ylabel(ax, 'Height - mm');
            title(ax, ['Flow Angle Field - ' num2str(sum(plot_logic)) ' Points'], 'FontWeight', 'bold');

        case 6

        density = field.rho;

            Density_Field = figure('Name', 'Density Field', 'NumberTitle', 'off', ...
            'Position', [100, 200, 1200, 320]);
            ax = axes('Parent', Density_Field);
            hold(ax, 'on');
            
            scatter(ax, X(plot_logic), Y(plot_logic), 20, density(plot_logic), 'filled');
            c = colorbar;               % Displays the color scale bar
            c.Label.String = 'Density'; % Labels the color key axis
            colormap('turbo'); 

            plot([0, (yMax - y_inj)/tand(beta + phi)], [y_inj, yMax], 'r--', 'LineWidth', 2); % Shock line
            plot([0, (yMax - y_inj)/tand(delta + phi)], [y_inj, yMax], 'g--', 'LineWidth', 2); % Slip line
            patch([x_infill_start, xMax, xMax, x_infill_start], [0, 0, y_inj, 0], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            patch([0, y_inj*tand(phi), 0, 0], [y_inj, 0, 0, y_inj], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle

            xlim(ax, [0, xMax]);
            ylim(ax, [0, chamber_length]);
            yticks(ax, 0:10:30);
            pbaspect(ax, [xMax, chamber_length + 10, 1]);
                
            grid(ax, 'on');
            box(ax, 'on');
                
            xlabel(ax, 'Axial Distance - mm (Unwrapped)');
            ylabel(ax, 'Height - mm');
            title(ax, ['Density Field - ' num2str(sum(plot_logic)) ' Points'], 'FontWeight', 'bold');

        case 7

        velocity = field.V;

            Velocity_Magnitude_Field = figure('Name', 'Velocity Field', 'NumberTitle', 'off', ...
            'Position', [100, 200, 1200, 320]);
            ax = axes('Parent', Velocity_Magnitude_Field);
            hold(ax, 'on');
            
            scatter(ax, X(plot_logic), Y(plot_logic), 20, velocity(plot_logic), 'filled');
            c = colorbar;               % Displays the color scale bar
            c.Label.String = 'Velocity Magnitude'; % Labels the color key axis
            colormap('turbo'); 

            plot([0, (yMax - y_inj)/tand(beta + phi)], [y_inj, yMax], 'r--', 'LineWidth', 2); % Shock line
            plot([0, (yMax - y_inj)/tand(delta + phi)], [y_inj, yMax], 'g--', 'LineWidth', 2); % Slip line
            patch([x_infill_start, xMax, xMax, x_infill_start], [0, 0, y_inj, 0], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            patch([0, y_inj*tand(phi), 0, 0], [y_inj, 0, 0, y_inj], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle

            xlim(ax, [0, xMax]);
            ylim(ax, [0, chamber_length]);
            yticks(ax, 0:10:30);
            pbaspect(ax, [xMax, chamber_length + 10, 1]);
                
            grid(ax, 'on');
            box(ax, 'on');
                
            xlabel(ax, 'Axial Distance - mm (Unwrapped)');
            ylabel(ax, 'Height - mm');
            title(ax, ['Velocity Magnitude Field - ' num2str(sum(plot_logic)) ' Points'], 'FontWeight', 'bold');

        otherwise
            warning('Unknown plot requested: %s', plots_requested(i));
    end
end

end
