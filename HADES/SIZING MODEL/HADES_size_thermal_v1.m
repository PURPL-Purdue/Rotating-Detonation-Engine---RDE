% Ansh Patel 
% RDE Chamber Hotfire Limit 
% 1D transient radial conduction in an outer wall
% with convection from hot detonation gas inside
% and convection to ambient air outside.
% Gas-side h is based on Ruan et al. (2026) 1D thermal model
% using recovery temperature and a Bartz-type Nusselt correlation.

% Flow annulus geometry (hot gas region)
ri_flow = 41.275e-3; % inner radius of annulus (1.625") [m]
ro_flow = 50.800e-3; % outer radius of annulus (2.000") [m]

% Solid wall geometry (what we are modeling conduction through)
wall_extra = 0.5 * 0.0254; % wall thickness [m] (0.5")
ri = ro_flow; % inner radius of solid wall [m]
ro = ro_flow + wall_extra;  % outer radius of solid wall [m]
L  = 0.093345;% axial length [m] (optional for heat per length)

T_init = 25;    % initial wall & ambient temperature [°C]
T_ambC = 25;    % ambient air outside [°C]
T_ambK = T_ambC + 273.15;   

% Inner-wall temperature limits to check (deg C)
T_limits_C = [700, 926.67, 1188.12];   % [°C]
nLimits    = numel(T_limits_C);
time_to_limit = NaN(nLimits,1);

% Design / target inner wall temperature
T_wallC_design = 1188.12; % [°C]
T_wallK_design = T_wallC_design + 273.15; % [K]

rho_wall = 7930; % kg/m^3
cp_wall = 500; % J/kg-K
k_wall = 15;     % W/m-K

%CEA DATA: 4 PHI VALUES, 1 CJ T EACH
% Equivalence ratios 
phi_vec = [0.8, 0.9, 1.0, 1.1];

% CJ temperatures corresponding to those phis (one per phi)
CJ_T_vec = [2790.67, 2914.30, 2990.66, 3019.33];  % [K]

% For each phi, we pick the CEA data at THAT CJ temperature column:
%   phi = 0.8 -> use column 1 of phi=0.8 block (T = 2790.67 K)
%   phi = 0.9 -> use column 2 of phi=0.9 block (T = 2914.30 K)
%   phi = 1.0 -> use column 3 of phi=1.0 block (T = 2990.66 K)
%   phi = 1.1 -> use column 4 of phi=1.1 block (T = 3019.33 K)

% ---- Static density rho_g [kg/m^3] at CJ for each phi ----
rho_g_vec = [ ...
    3.3327, ...  % phi = 0.8 @ 2790.67 K
    3.1250, ...  % phi = 0.9 @ 2914.30 K
    2.9747, ...  % phi = 1.0 @ 2990.66 K
    2.8724  ...  % phi = 1.1 @ 3019.33 K
];

% Cp with equilibrium reactions [kJ/(kg·K)] at CJ, then to J/(kg·K)
cp_g_kJ_vec = [ ...
    2.3336, ...  % phi = 0.8, col1
    2.7540, ...  % phi = 0.9, col2
    3.1354, ...  % phi = 1.0, col3
    3.2067  ...  % phi = 1.1, col4
];

cp_g_vec = 1000 .* cp_g_kJ_vec;   % J/(kg·K)

% Gamma (specific heat ratio) at CJ 
gamma_g_vec = [ ...
    1.1924, ...  % phi = 0.8, col1
    1.1780, ...  % phi = 0.9, col2
    1.1693, ...  % phi = 1.0, col3
    1.1700  ...  % phi = 1.1, col4
];

% Prandtl number WITH equilibrium reactions at CJ 
Pr_t_g_vec = [ ...
    0.5242, ...  % phi = 0.8, col1
    0.4850, ...  % phi = 0.9, col2
    0.4988, ...  % phi = 1.0, col3
    0.5123  ...  % phi = 1.1, col4
];

%Viscosity at CJ [millipoise] - Pa·s 
visc_mP_vec = [ ...
    0.89074, ... % phi = 0.8, col1
    0.91898, ... % phi = 0.9, col2
    0.93556, ... % phi = 1.0, col3
    0.94054  ... % phi = 1.1, col4
];
point8_CJ = [0.6, 0.57, 0.57, 0.60,];
mu_g_vec = visc_mP_vec .* 1e-4;   % Pa·s (1 mP = 1e-4 Pa·s)


% phiIndex = 1 -> 0.8
% phiIndex = 2 -> 0.9
% phiIndex = 3 -> 1.0
% phiIndex = 4 -> 1.1
phiIndex = 2;             
phi = phi_vec(phiIndex);

% CJ temperature for this phi
T_CJ = CJ_T_vec(phiIndex); % [K]

% Grab gas properties for this phi at its CJ temperature
rho_g   = rho_g_vec(phiIndex);
cp_g = cp_g_vec(phiIndex);
gamma_g = gamma_g_vec(phiIndex);
Pr_t_g  = Pr_t_g_vec(phiIndex);
mu_g    = mu_g_vec(phiIndex);

T_gasK = point8_CJ(phiIndex) * T_CJ; % [K]

fprintf('\n=== Gas state (phi = %.2f) using its CJ T ===\n', phi);
fprintf('CJ T(phi) = %.2f K\n', T_CJ);
fprintf('T_gasK (0.6CJ) = %.2f K (%.2f °C)\n', T_gasK, T_gasK - 273.15);
fprintf('rho_g  = %.4f kg/m^3\n', rho_g);
fprintf('cp_g   = %.2f J/(kg·K)\n', cp_g);
fprintf('gamma  = %.4f\n', gamma_g);
fprintf('Pr_t,g = %.4f\n', Pr_t_g);
fprintf('mu_g   = %.4e Pa·s\n', mu_g);


%FILM TEMPERATURE & GAS-SIDE h_inner (Ruan/Bartz)
% Film temperature using inner wall temperature
T_filmK = 0.5 * (T_gasK + T_wallK_design); % [K]

% Adiabatic-mean density at film T (simple ideal-gas scaling)
rho_am_g = rho_g * (T_gasK / T_filmK); % kg/m^3


mu_am_g  = mu_g; % Pa·s

% Axial gas velocity & Mach (can refine later if you want them phi-dependent)
U_g = 1952.6 * 0.6; % m/s, bulk gas velocity (placeholder)
M_g = 0.3; % Mach number at annulus (placeholder)

% Recovery temperature (Ruan Eq. (3),(4))
theta_rec = Pr_t_g^(1/3); % theta = Pr^(1/3)
T_rec_g = T_gasK * (1 + theta_rec * (gamma_g - 1)/2 * M_g^2); % [K]

% Thermal conductivity via μ cp / Pr
k_gas = mu_g * cp_g / Pr_t_g; % W/(m·K)

% Hydraulic diameter for annulus [m]
D_h = 2 * (ro_flow - ri_flow);

% Reynolds number
Re_g = rho_g * U_g * D_h / mu_g;

% Bartz-style Nusselt (Ruan Eq. (5))
Nu_g = 0.026 * Re_g^0.8 * Pr_t_g^0.4 * (rho_am_g / rho_g)^0.8 * (mu_am_g / mu_g)^0.2;

% Gas-side convective h on inner wall
h_inner = Nu_g * k_gas / D_h;   % W/m^2-K

% Outer convection coefficient (ambient air)
h_outer = 15;  % W/m^2-K (approximate)

fprintf('T_wall_design = %.1f K (%.1f °C)\n', T_wallK_design, T_wallC_design);
fprintf('T_filmK = %.1f K\n', T_filmK);
fprintf('T_rec_g  = %.1f K (%.1f °C)\n', T_rec_g, T_rec_g - 273.15);
fprintf('Re_g = %.3e, Nu_g = %.1f, h_inner = %.1f W/m^2-K\n', Re_g, Nu_g, h_inner);

%  7. PACKAGE PARAMETERS FOR PDE
params.rho = rho_wall;
params.cp = cp_wall;
params.k = k_wall;
params.h_inner = h_inner;
params.h_outer = h_outer;

% Use gas recovery temp (in °C) as the effective gas-side boundary temp
params.T_gas = T_rec_g - 273.15;  
params.T_amb = T_ambC;
params.T_init = T_init;

%SPATIAL / TIME DISCRETIZATION (PDEPE)
Nr = 120;
r  = linspace(ri, ro, Nr); % radial grid [m]

tEnd = 2; % total simulated time [s]
Nt  = 400;
t = linspace(0, tEnd, Nt); % time grid

m = 1;  % cylindrical symmetry: radial coordinate r
% Equation: ρcp dT/dt = (1/r) d/dr( k r dT/dr )   in ri <= r <= ro
% pdepe form: c*u_t = 1/x^m d/dx( x^m f ) + s, with m = 1, f = k*dT/dr, s = 0
 
sol = pdepe(m, @(r,t,u,dudr) heatPDE(r,t,u,dudr,params), @(r)icFun(r,params),@(rl,ul,rr,ur,t) bcFun(rl,ul,rr,ur,t,params),r, t);
T = sol(:,:,1);   % T(t_index, r_index) in °C

%INNER WALL TEMPERATURE VS TIME
T_inner = T(:,1);% [°C] at r = ri
T_inner_F = T_inner * 9/5 + 32;  % [°F]
T_limits_F = T_limits_C * 9/5 + 32;  % [°F] thresholds

% vRUNTIME UNTIL EACH LIMIT IS REACHED (°C)
for i = 1:nLimits
    Tcrit = T_limits_C(i);

    if (T_inner(1) >= Tcrit) || all(T_inner < Tcrit)
        time_to_limit(i) = NaN;  % never reaches or already above at t=0
    else
        time_to_limit(i) = interp1(T_inner, t, Tcrit, 'linear');
    end
end


fprintf('\n=== Max runtime based on INNER WALL temperature (phi = %.2f) ===\n', phi);
for i = 1:nLimits
    if isnan(time_to_limit(i))
        fprintf('Inner wall never reaches %.2f °C within %.2f s simulated.\n', ...
            T_limits_C(i), tEnd);
    else
        fprintf('Time for inner wall to reach %.2f °C (%.2f °F): %.5f s\n', ...
            T_limits_C(i), T_limits_F(i), time_to_limit(i));
    end
end



% 1) Inner wall temperature vs time (°F)
figure;
plot(t, T_inner_F, 'Color', '#9100FF', 'LineWidth', 5);
hold on;
for i = 1:nLimits
    yline(T_limits_F(i), '--', sprintf('%.2f °F', T_limits_F(i)));
end
xlabel('Time [s]');
ylabel('Inner Wall Temperature [°F]');
title(sprintf('Inner Wall Temperature vs Time (\\phi = %.2f)', phi));
grid on;

% 2) Radial temperature profiles at a few times (°F)
times_to_plot = [0.01, 0.05, 0.1, 0.2, 0.5]; % s 
figure; hold on;
for tp = times_to_plot
    [~, idx] = min(abs(t - tp));
    T_profile_F = T(idx,:) * 9/5 + 32;   % °C -> °F
    plot(r*1000, T_profile_F, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('t = %.3f s', t(idx)));
end
xlabel('Radius [mm]');
ylabel('Temperature [°F]');
title(sprintf('Radial temperature profiles in wall (\\phi = %.2f)', phi));
legend('show', 'Location', 'best');
grid on;

% LOCAL FUNCTIONS FOR PDEPE 

function [c,f,s] = heatPDE(~,~,~,dudr,params)
    % ρcp dT/dt = (1/r) d/dr (k r dT/dr)
    c = params.rho * params.cp;
    f = params.k * dudr;
    s = 0;
end

function u0 = icFun(~,params)
    % Initial condition: uniform T_init
    u0 = params.T_init;
end

function [pl,ql,pr,qr] = bcFun(~,ul,~,ur,~,params)
    % Inner boundary r = ri:
    %   -k dT/dr = h_inner (T_gas - T_wall)
    pl = params.h_inner * (params.T_gas - ul);
    ql = 1;

    % Outer boundary r = ro:
    %   -k dT/dr = h_outer (T_wall - T_amb)
    pr = params.h_outer * (ur - params.T_amb);
    qr = 1;
end
