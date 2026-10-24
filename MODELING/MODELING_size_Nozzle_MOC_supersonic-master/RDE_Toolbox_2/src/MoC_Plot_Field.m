function [] = MoC_Plot_Field(field, geometry, plots_requested, plot_options)
% MoC_Plot_Field Plots the unwrapped Method of Characteristics net and RDE flow domain.
%
% Inputs:
%   field: Structure containing the solved characteristic points and properties
%   geometry: Structure containing domain parameters (xMax, yMax, y_inj, x_inj, phi, etc.)
%   plots_requested: Numeric vector specifying which plots to generate
%       1 = Method of Characteristics Net
%       2 = Static Temperature Field [K]
%       3 = Static Pressure Field [bar]
%       4 = Mach Number Field
%       5 = Flow Angle (Theta) Field [deg]
%       6 = Density Field [kg/m^3]
%       7 = Velocity Magnitude Field [m/s]
%       8 = Characteristic Lines (products: C+ red, C- blue; slip-shock: green)
%           Dimensionless x/L and y/L, L = chamber circumference (xMax); no grid.
%   Uses Parent_C_Plus/Minus connectivity and X_unwrapped marching coordinates.

%   plot_options: Vector of numbers specifying which plot options to include
%       1 = Show Slip Line
%       2 = Show Shock Line
%       3 = Show Infill Region
%       4 = Show Detonation Triangle

%% Extract parameters

Y = field.Y;
X = field.X;
Region = field.Region;
plot_logic = field.plot_logic; % Logical index for points within the chamber (Region 1, 2, and 3)
datapoint_indices = field.datapoint_indices; % Index numbers for all points in the field
y_inj = geometry.y_inj;
phi = geometry.phi;
beta = geometry.beta;
delta = geometry.delta;
yMax = geometry.yMax;
xMax = geometry.xMax;
chamber_length = geometry.chamber_length;
x_infill_start = geometry.x_infill_start;

slipline = any(plot_options == 1);
shockline = any(plot_options == 2);
infillregion = any(plot_options == 3);
detonationtriangle = any(plot_options == 4);

%% Plot 

for i = 1:length(plots_requested)
    switch plots_requested(i)
        case {1, 8}

            net_title = 'Method of Characteristics Net';
            length_scale = 1;
            if plots_requested(i) == 8
                net_title = 'Method of Characteristics Lines';
                length_scale = xMax;
            end
            base_figure = figure('Name', net_title, 'NumberTitle', 'off', ...
            'Position', [100, 200, 1200, 320]);
            ax = axes('Parent', base_figure);
            hold(ax, 'on');

            idx = find(plot_logic);
            if plots_requested(i) == 8
                plotCharacteristicLines(ax, field, xMax);
            else
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
            end

            plot(ax, [0, (yMax - y_inj)/tand(beta + phi)]/length_scale, [y_inj, yMax]/length_scale, 'r--', 'LineWidth', 2); % Shock line
            plot(ax, [0, (yMax - y_inj)/tand(delta + phi)]/length_scale, [y_inj, yMax]/length_scale, 'g--', 'LineWidth', 2); % Slip line
            patch(ax, [x_infill_start, xMax, xMax, x_infill_start]/length_scale, [0, 0, y_inj, 0]/length_scale, 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            patch(ax, [0, y_inj*tand(phi), 0, 0]/length_scale, [y_inj, 0, 0, y_inj]/length_scale, 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle

            xlim(ax, [0, xMax]/length_scale);
            ylim(ax, [0, chamber_length]/length_scale);
            daspect(ax, [1 1 1]);
            box(ax, 'off');
            if plots_requested(i) == 8
                grid(ax, 'off');
                xlabel(ax, 'x/L (L = chamber circumference)');
                ylabel(ax, 'y/L');
            else
                yticks(ax, 0:10:30);
                grid(ax, 'on');
                xlabel(ax, 'Axial Distance - mm (Unwrapped)');
                ylabel(ax, 'Height - mm');
            end
            title(ax, [net_title ' - ' num2str(length(X)) ' Points'], 'FontWeight', 'bold');
            
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

            if shockline
                plot([0, (yMax - y_inj)/tand(beta + phi)], [y_inj, yMax], 'r--', 'LineWidth', 2); % Shock line
            end
            if slipline
                plot([0, (yMax - y_inj)/tand(delta + phi)], [y_inj, yMax], 'g--', 'LineWidth', 2); % Slip line
            end
            if infillregion
                patch([x_infill_start, xMax, xMax, x_infill_start], [0, 0, y_inj, 0], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            end
            if detonationtriangle
                patch([0, y_inj*tand(phi), 0, 0], [y_inj, 0, 0, y_inj], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle
            end

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

            if shockline
                plot([0, (yMax - y_inj)/tand(beta + phi)], [y_inj, yMax], 'r--', 'LineWidth', 2); % Shock line
            end
            if slipline
                plot([0, (yMax - y_inj)/tand(delta + phi)], [y_inj, yMax], 'g--', 'LineWidth', 2); % Slip line
            end
            if infillregion
                patch([x_infill_start, xMax, xMax, x_infill_start], [0, 0, y_inj, 0], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            end
            if detonationtriangle
                patch([0, y_inj*tand(phi), 0, 0], [y_inj, 0, 0, y_inj], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle
            end

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

            if shockline
                plot([0, (yMax - y_inj)/tand(beta + phi)], [y_inj, yMax], 'r--', 'LineWidth', 2); % Shock line
            end
            if slipline
                plot([0, (yMax - y_inj)/tand(delta + phi)], [y_inj, yMax], 'g--', 'LineWidth', 2); % Slip line
            end
            if infillregion
                patch([x_infill_start, xMax, xMax, x_infill_start], [0, 0, y_inj, 0], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            end
            if detonationtriangle
                patch([0, y_inj*tand(phi), 0, 0], [y_inj, 0, 0, y_inj], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle
            end

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

            if shockline
                plot([0, (yMax - y_inj)/tand(beta + phi)], [y_inj, yMax], 'r--', 'LineWidth', 2); % Shock line
            end
            if slipline
                plot([0, (yMax - y_inj)/tand(delta + phi)], [y_inj, yMax], 'g--', 'LineWidth', 2); % Slip line
            end
            if infillregion
                patch([x_infill_start, xMax, xMax, x_infill_start], [0, 0, y_inj, 0], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            end
            if detonationtriangle
                patch([0, y_inj*tand(phi), 0, 0], [y_inj, 0, 0, y_inj], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle
            end

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

            if shockline
                plot([0, (yMax - y_inj)/tand(beta + phi)], [y_inj, yMax], 'r--', 'LineWidth', 2); % Shock line
            end
            if slipline
                plot([0, (yMax - y_inj)/tand(delta + phi)], [y_inj, yMax], 'g--', 'LineWidth', 2); % Slip line
            end
            if infillregion
                patch([x_infill_start, xMax, xMax, x_infill_start], [0, 0, y_inj, 0], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            end
            if detonationtriangle
                patch([0, y_inj*tand(phi), 0, 0], [y_inj, 0, 0, y_inj], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle
            end

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

            if shockline
                plot([0, (yMax - y_inj)/tand(beta + phi)], [y_inj, yMax], 'r--', 'LineWidth', 2); % Shock line
            end
            if slipline
                plot([0, (yMax - y_inj)/tand(delta + phi)], [y_inj, yMax], 'g--', 'LineWidth', 2); % Slip line
            end
            if infillregion
                patch([x_infill_start, xMax, xMax, x_infill_start], [0, 0, y_inj, 0], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Infill triangle
            end
            if detonationtriangle
                patch([0, y_inj*tand(phi), 0, 0], [y_inj, 0, 0, y_inj], 'c', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % Detonation triangle
            end

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

% Sum bs slop to accurately draw the char lines when plot 8 is requested. 
function plotCharacteristicLines(ax, field, xMax)
% Draw actual parent-to-child segments, batching each color into one object.
% Marching coordinates retain which side of the periodic seam each point is on.
if isfield(field, 'X_unwrapped')
    x = field.X_unwrapped(:);
else
    x = field.X(:); % Also accepts a field before its coordinates are wrapped.
end
y = field.Y(:);
region = field.Region(:);
visible = field.plot_logic(:);
line_x = {[], [], []};
line_y = {[], [], []};
families = {'Parent_C_Plus', 'Parent_C_Minus'};
for family = 1:2
    parents = field.(families{family})(:);
    children = find(isfinite(parents) & parents >= 1 & ...
        parents <= numel(x) & parents == fix(parents));
    parents = parents(children);
    keep = (visible(parents) | visible(children)) & ...
        isfinite(x(parents)) & isfinite(x(children)) & ...
        isfinite(y(parents)) & isfinite(y(children));
    parents = parents(keep);
    children = children(keep);
    % An incident product ray ending at a shock seed is still a product ray.
    % Rays emitted by the shock seeds and the slip-shock net are green.
    slip_shock = region(parents) == 3;
    for in_slip_shock = 0:1
        select = slip_shock == in_slip_shock;
        color_index = family;
        if in_slip_shock
            color_index = 3;
        end
        [sx, sy] = wrapCharacteristicSegments(x(parents(select)), ...
            y(parents(select)), x(children(select)), y(children(select)), xMax);
        line_x{color_index} = [line_x{color_index}; sx];
        line_y{color_index} = [line_y{color_index}; sy];
    end
end
colors = {'r', 'b', 'g'};
labels = {'Product C+', 'Product C-', 'Slip-shock characteristics'};
handles = gobjects(1, 3);
for k = 1:3
    if isempty(line_x{k})
        line_x{k} = NaN;
        line_y{k} = NaN;
    end
    % Normalize both directions by the same length to retain the geometry.
    handles(k) = plot(ax, line_x{k}/xMax, line_y{k}/xMax, 'Color', colors{k}, ...
        'LineWidth', 0.75, 'DisplayName', labels{k});
end
legend(ax, handles, labels, 'Location', 'best', 'AutoUpdate', 'off');
end

function [sx, sy] = wrapCharacteristicSegments(x1, y1, x2, y2, period)
% Clip each segment into periodic strips before mapping it into [0, period].
% NaNs separate segments, preventing connections between unrelated rays.
sx = zeros(0, 1);
sy = zeros(0, 1);
if isempty(x1)
    return
end
dx = x2 - x1;
dy = y2 - y1;
for strip = floor(min([x1; x2])/period):floor(max([x1; x2])/period)
    left = strip*period;
    right = left + period;
    t0 = zeros(size(x1));
    t1 = ones(size(x1));
    moving = dx ~= 0;
    a = (left-x1(moving))./dx(moving);
    b = (right-x1(moving))./dx(moving);
    t0(moving) = max(0, min(a, b));
    t1(moving) = min(1, max(a, b));
    keep = t0 < t1 & max(x1, x2) >= left & min(x1, x2) <= right;
    xs = [x1(keep)+t0(keep).*dx(keep)-left, ...
          x1(keep)+t1(keep).*dx(keep)-left, nan(nnz(keep), 1)]';
    ys = [y1(keep)+t0(keep).*dy(keep), ...
          y1(keep)+t1(keep).*dy(keep), nan(nnz(keep), 1)]';
    sx = [sx; xs(:)];
    sy = [sy; ys(:)];
end
end
