clc; close;

d = 0:0.5:4; % depth (m)
F = 24525 .* sqrt((d .^ 4) + (d .^ 6) ./ 36); %hydrostatic force (N)
plot(h, F, "LineWidth", 1.5);
xlabel("depth (m)");
ylabel("force (N)");
title("Hydrostatic Force vs Depth");
grid on;
