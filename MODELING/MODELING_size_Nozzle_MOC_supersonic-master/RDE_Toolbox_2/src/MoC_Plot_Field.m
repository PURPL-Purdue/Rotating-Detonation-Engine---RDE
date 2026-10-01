function [fieldFigure, netFigure] = MoC_Plot_Field(field, geometry, xMax, yMax)
% Plots the unwrapped field and characteristic net.

% Inputs:
%   field: structure containing the solved points and their properties
%   geometry: structure containing the chamber dimensions and other geometric parameters
%   xMax, yMax: maximum dimensions of the plotting domain

% Outputs:
%   fieldFigure: handle to the figure displaying the flow field
%   netFigure: handle to the figure displaying the characteristic net

% Taking chamber dimensions and plotting graph domain: 

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



