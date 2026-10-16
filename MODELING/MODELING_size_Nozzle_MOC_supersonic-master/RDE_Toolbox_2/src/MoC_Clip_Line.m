function [xc, yc, visible] = MoC_Clip_Line(x1, y1, x2, y2, xMax, yMax, y_inj, beta)
% MoC_Clip_Line Clips a characteristic line segment to the RDE domain boundaries.
%
% Boundaries for Cycle 1:
%   - Left boundary: x >= 0
%   - Bottom boundary: y >= 0
%   - Top boundary: y <= yMax (chamber length exit)
%   - Downstream boundary:
%       * For y <= y_inj: vertical detonation / infill rear boundary (x <= xMax)
%       * For y > y_inj:  oblique shock boundary (x <= xMax + (y - y_inj)/tand(beta))

    %% 1. Basic Box Clipping (x >= 0, 0 <= y <= yMax)
    dx = x2 - x1;
    dy = y2 - y1;
    t0 = 0;
    t1 = 1;

    % Clip to y >= 0
    if dy == 0
        if y1 < 0, xc = []; yc = []; visible = false; return; end
    else
        t = -y1 / dy;
        if dy > 0, t0 = max(t0, t); else, t1 = min(t1, t); end
    end

    % Clip to y <= yMax (top boundary)
    if dy == 0
        if y1 > yMax, xc = []; yc = []; visible = false; return; end
    else
        t = (yMax - y1) / dy;
        if dy > 0, t1 = min(t1, t); else, t0 = max(t0, t); end
    end

    % Clip to x >= 0
    if dx == 0
        if x1 < 0, xc = []; yc = []; visible = false; return; end
    else
        t = -x1 / dx;
        if dx > 0, t0 = max(t0, t); else, t1 = min(t1, t); end
    end

    % If box clipping leaves nothing
    if t0 >= t1
        xc = []; yc = []; visible = false; return;
    end

    % If no shock boundary is provided, perform standard box clipping at x <= xMax
    if nargin < 7 || isempty(y_inj) || isempty(beta)
        if dx == 0
            if x1 > xMax, xc = []; yc = []; visible = false; return; end
        else
            t = (xMax - x1) / dx;
            if dx > 0, t1 = min(t1, t); else, t0 = max(t0, t); end
        end
        if t0 >= t1
            xc = []; yc = []; visible = false; return;
        end
        xc = [x1 + t0 * dx, x1 + t1 * dx];
        yc = [y1 + t0 * dy, y1 + t1 * dy];
        visible = true;
        return;
    end

    %% 2. Downstream Shock / Detonation Boundary Clipping
    xa = x1 + t0 * dx; ya = y1 + t0 * dy;
    xb = x1 + t1 * dx; yb = y1 + t1 * dy;

    tb = tand(beta);
    % Downstream boundary definition:
    % x_bound(y) = xMax for y <= y_inj
    % x_bound(y) = xMax + (y - y_inj)/tb for y > y_inj
    x_bnd = @(y) xMax + max(0, y - y_inj) / tb;

    in_a = (xa <= x_bnd(ya) + 1e-6);
    in_b = (xb <= x_bnd(yb) + 1e-6);

    if in_a && in_b
        xc = [xa, xb];
        yc = [ya, yb];
        visible = true;
        return;
    elseif ~in_a && ~in_b
        xc = [];
        yc = [];
        visible = false;
        return;
    end

    % Exactly one endpoint is inside and one is outside.
    % Find parameter t along segment [xa, ya] -> [xb, yb] where it crosses downstream boundary.
    t_int = [];
    del_x = xb - xa;
    del_y = yb - ya;

    % 1. Try vertical boundary segment: x = xMax (for y in [0, y_inj])
    if abs(del_x) > 1e-9
        tv = (xMax - xa) / del_x;
        if tv >= -1e-6 && tv <= 1 + 1e-6
            yv = ya + tv * del_y;
            if yv >= -1e-6 && yv <= y_inj + 1e-6
                t_int = max(0, min(1, tv));
            end
        end
    end

    % 2. Try oblique shock boundary: f(x,y) = tb*(x - xMax) - (y - y_inj) = 0 (for y >= y_inj)
    if isempty(t_int)
        fa = tb * (xa - xMax) - (ya - y_inj);
        fb = tb * (xb - xMax) - (yb - y_inj);
        denom = fb - fa;
        if abs(denom) > 1e-9
            ts = -fa / denom;
            if ts >= -1e-6 && ts <= 1 + 1e-6
                ys = ya + ts * del_y;
                if ys >= y_inj - 1e-6 && ys <= yMax + 1e-6
                    t_int = max(0, min(1, ts));
                end
            end
        end
    end

    if isempty(t_int)
        % Fallback for extreme corner precision
        if in_a
            xc = [xa, xa]; yc = [ya, ya];
        else
            xc = [xb, xb]; yc = [yb, yb];
        end
        visible = true;
        return;
    end

    x_cross = xa + t_int * del_x;
    y_cross = ya + t_int * del_y;

    if in_a
        xc = [xa, x_cross];
        yc = [ya, y_cross];
    else
        xc = [x_cross, xb];
        yc = [y_cross, yb];
    end
    visible = true;
end

