%% The big FFC

%{
NOTES
    - Requires Optimization toolbox add-on
    - Oxygen is first column, methane second
    - Don't guess 1 for M_G,  calc needs either subsonic or supersonic guess
    - All other guesses used only for Darcy friction calculations,
    iterative convergence with Fanno fixes output values
    - Loops utilized for all solvers because they don't work well with vectors
%}

clear; clc;

%% ---------- General Inputs ----------

Thrust = 1780;          % (N)
Isp = 233.37;           % Specific impulse, from CEA (s)
Pc = 1e6;               % Chamber pressure (Pa)
cstar = 1806.2;         % Characteristic velocity, from CEA (m/s)
phi = 1.3;              % Equivalence ratio 
Cd = 0.6;               % Discharge coefficient
elements = 36;          % Number of injector elements
Imp_angle = 45;         % Injector impingement angle, off axial (deg)
ope = [1 2];            % Orifices per injector element
inj_OD = [2.5 1];       % Injector diameter (mm)
T = 283.15;             % Satic inlet temperature
Rmid = 22.3524953;      % Average annulus diameter
M_G = [0.5 0.5];        % Speed regime guess for FFC, used in Re calc as well
P_G = [3e6 3e6];        % Pressure guess for Re (Pa)
T_stagG = [300 300];    % Temperature guess for Re (K)
f_G = [0.02 0.02];      % Darcy friction factor guess


% Fanno Flow inputs

%L = [9.48944 9.56818];  % Tube length 
L = [9.5 9.5];
D = inj_OD;             % Tube diameter 

% Constants

g = 9.81;           % (m/s^2)
gamma = [1.4 1.31]; % specific heat ratio
R = [259.8 518.3];  % Gas constants (J/kg-K)

mu_ref = [2.07e-5 1.1e-5];        % Reference viscosity @ reference temp (Pa-s)
T_ref = [291 291];                % Reference termperature (K)
S_ref = [127 148];                % Reference Sutherland Constant (K)
epsilon = 0.015;                  % Stainless Steel surface roughness (mm)


%% ---------- calc ----------

OF_ratio = 4./phi;                                  % Oxygen to fuel ratio
Y = [(OF_ratio/(1 + OF_ratio)) (1/(1 + OF_ratio))]; % Propellant mass fraction
m_dot = (Y .* Thrust) ./ (Isp .* g);                % Propellant Mass flow rate (kg/s)
m_dotT = Thrust ./ (Isp .* g);                      % Total mass flow rate (kg/s)
A_inj = (inj_OD./2000).^2 .* pi .* elements .* ope; % Injector area (m^2)
A_t = ((m_dotT .* cstar) ./ Pc) * 1e6;              % Throat area (mm^2)
A_gap = A_t ./ (Rmid .* 2*pi);                      % Annulus gap (mm)
Rinner = Rmid - A_gap ./ 2;                         % Inner throat radius (mm)
Router = Rmid + A_gap ./ 2;                         % Outer throat radius (mm)


% ---------- Iteration Setup ----------

tol = 1e-6;
max_iter = 100;
err = 1;
iter = 0;

P2 = P_G;        % initial guess (Pa)
T_stag = T_stagG;        % initial guess (K)

while err > tol && iter < max_iter
    
    iter = iter + 1;
    
    % ---- Update viscosity (Sutherland) ----
    mu = mu_ref .* (T_ref + S_ref) ./ (T_stag + S_ref) ...
         .* (T_stag ./ T_ref).^(3/2);
    
    % ---- Reynolds number ----
    Re = (P2 .* M_G .* sqrt(gamma .* R .* T) .* L) ...
         ./ (R .* T .* mu);
    
  

    % ---- Loop over propellants ----
    for i = 1:2
    
    % ---- Friction factor via Colebrook ----
        colebrook = @(ff) 1./sqrt(ff) + ...
            2*log10( epsilon/(3.7*inj_OD(i)) + 2.51/(Re(i)*sqrt(ff)) );
    
        f(i) = fsolve(colebrook, f_G(i), optimoptions('fsolve','Display','off'));
    
    % ---- Solve Fanno for inlet Mach ----
        fanno_fun = @(M) ...
            (1 - M^2)/(gamma(i)*M^2) + ((gamma(i)+1)/(2*gamma(i))) ...
            * log(((gamma(i)+1)*M^2)/(2+(gamma(i)-1)*M^2)) - f(i)*L(i)/D(i);
    
        M_TubeIN(i) = fsolve(fanno_fun, M_G(i), optimoptions('fsolve','Display','off'));
    
        if M_TubeIN(i) <= 0
            error('Solver returned nonphysical Mach number.');
        end
    
    % ---- Solve for upstream pressure P01 using fzero ----
        P01_fun = @(P0) Cd*A_inj(i)*P0*sqrt((gamma(i)/(R(i)*T))*(P2(i)/P0)^(2/gamma(i)) ...
                *(1-(P2(i)/P0)^((gamma(i)-1)/gamma(i)))) - m_dot(i);

        a(i) = P2(i)*1.01;
        b(i) = 1e7*a(i);  % Fixed upper bound
        
        if P01_fun(a(i))*P01_fun(b(i)) < 0
            P01(i) = fzero(P01_fun, [a(i) b(i)]);
        else
            P01(i) = NaN;
        end
    end
    
    % ---- Update stagnation properties ----
    T_stag_new = T .* (1 + ((gamma-1)/2) .* M_TubeIN.^2);
    
   
    P_new = Pc .*(1./M_TubeIN).*(1./sqrt((2./(gamma+1))...
        .*(1+((gamma-1)./2).*M_TubeIN.^2)));
    

    % ---- Convergence check (pressure only) ----
    err = max(abs((P_new - P2)./P2));
    
    % Update guesses
    P2 = P_new;
    T_stag = T_stag_new;
    
end

if iter == max_iter
    warning('Iteration did not fully converge');
end

%% P01 solver convergence diagnostic

%{
x = linspace(P2(1), 5*P2(2), 1e3);
y1 = Cd .* A_inj(1) .* x .* sqrt((gamma(1) ./ (R(1) .* T)) .* (P2(1) ./ x).^(2 ./ gamma(1)) ...
        .* (1 - (P2(1) ./ x).^((gamma(1) - 1) ./ gamma(1))));
y2 = Cd .* A_inj(2) .* x .* sqrt((gamma(2) ./ (R(2) .* T)) .* (P2(2) ./ x).^(2 ./ gamma(2)) ...
        .* (1 - (P2(2) ./ x).^((gamma(2) - 1) ./ gamma(2))));
y12 = linspace(m_dot(1), m_dot(1), 1e3);

figure
plot(x, y1);
hold on
plot(x, y12);
plot(x, y2);
%}

%% ---------- Output ----------

fprintf('\n');
fprintf('====================================================\n');
fprintf('              FLOW CALCULATION SUMMARY              \n');
fprintf('====================================================\n\n');

fprintf('--- Mass Flow Rates ---\n');
fprintf('Total Mass Flow Rate            : %.3f kg/s\n', m_dotT);
fprintf('Oxygen Mass Flow Rate           : %.3f kg/s\n', m_dot(1));
fprintf('Methane Mass Flow Rate          : %.3f kg/s\n\n', m_dot(2));

fprintf('--- Annular Throat Geometry ---\n');
fprintf('Outer Radius                    : %.2f mm\n', Router);
fprintf('Inner Radius                    : %.2f mm\n', Rinner);
fprintf('Annulus Gap                     : %.2f mm\n\n', A_gap);

fprintf('--- Fanno Output & Pressure Correlations ---\n');
fprintf('Oxygen Inlet Mach               : %.2f \n', M_TubeIN(1));
fprintf('Methane Inlet Mach              : %.2f \n', M_TubeIN(2));
fprintf('P2 Oxygen Downstream Pressure   : %.2f bar | %.2f psi\n', P2(1)/1e5, (P2(1)*14.504)/1e5);
fprintf('P2 Methane Downstream Pressure  : %.2f bar | %.2f psi\n', P2(2)/1e5, (P2(2)*14.504)/1e5);
fprintf('P01 Oxygen Upstream Pressure    : %.2f bar | %.2f psi\n' , P01(1)/1e5, (P01(1)*14.504)/1e5);
fprintf('P01 Methane Upstream Pressure   : %.2f bar | %.2f psi\n\n' , P01(2)/1e5, (P01(2)*14.504)/1e5);

fprintf('====================================================\n\n');



