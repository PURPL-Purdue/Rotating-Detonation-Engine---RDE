function ax = Setup_Graph_Domain(xMax, yMax)
%SETUP_GRAPH_DOMAIN Create the empty plotting domain for the RDE graph.


%
%   Inputs:
%       xMax - domain extent along x (same units your computation uses)
%       yMax - domain extent along y
%
%   Output:
%       ax   - axes handle to plot onto in later post-processing steps

fig = figure('Name', 'RDE Flow Field', 'NumberTitle', 'off');
ax = axes('Parent', fig);

xlim(ax, [0, xMax]);
ylim(ax, [0, yMax]);

% Lock the plot box's shape to match the data ratio (xMax:yMax) so
% one unit in x and one unit in y are the same length on screen.
% This matters for MoC: characteristic/wave angles only look correct
% under true 1:1 scaling. Using pbaspect (not axis equal) gets the
% same 1:1 scaling WITHOUT MATLAB expanding xlim/ylim to solve for it.
pbaspect(ax, [xMax, yMax, 1]);

xlabel(ax, 'x (mm)');      % TODO: set real units once decided (m, mm, deg...)
ylabel(ax, 'y (mm)');
title(ax, 'RDE Domain');

grid(ax, 'on');
hold(ax, 'on');       % so later functions can add to this axes
% instead of overwriting it
end