%% The big FFC

%{
NOTES
    - Requires Optimization toolbox add-on
    - Oxygen is first column, methane second
    - Don't guess 1 for M_G,  calc needs either subsonic or supersonic guess
    - All other guesses used only for Darcy friction calculations,
    iterative convergence with Fanno fixes output values
    - P01 Solver Diagnositic is toggled above the output
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

P2 = P_G;           % initial guess (Pa)
T_stag = T_stagG;   % initial guess (K)

while err > tol && iter < max_iter
    
    iter = iter + 1;
    
    % ---- Update viscosity (Sutherland) ----
    mu = mu_ref .* (T_ref + S_ref) ./ (T_stag + S_ref) ...
         .* (T_stag ./ T_ref).^(3/2);
    
    % ---- Reynolds number ----
    Re = (P2 .* M_G .* sqrt(gamma .* R .* T) .* (D./1000)) ...
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
            error('ERROR: Fanno Solver returned nonphysical Mach number.');
        end
    
    % ---- Solve for upstream pressure P01 using fzero ----
        P01_fun = @(P0) Cd*A_inj(i)*P0*sqrt((2*gamma(i)/(R(i)*T*(gamma(i)-1)))*(P2(i)/P0)^(2/gamma(i)) ...
                *(1-(P2(i)/P0)^((gamma(i)-1)/gamma(i)))) - m_dot(i);

        a(i) = P2(i)*1.01;
        b(i) = a(i) + 1e8;  % 1000 bar (1e8 Pa) Solution Window
        
        if P01_fun(a(i))*P01_fun(b(i)) < 0
            P01(i) = fzero(P01_fun, [a(i) b(i)]);
        else
            error('ERROR: P01 Solver Solution Window Does Not Include A Valid Root');
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

%% Toggle P01 Solver Diagnostic

%{

for i = 1:2
    x(i,:) = P2(i) * (1 + logspace(-6, 0.25, 3e3));
    ym(i,:) = linspace(m_dot(i), m_dot(i), 3e3);

    y(i,:) = Cd.*A_inj(i).*x(i,:).*sqrt((2.*gamma(i)./(R(i).*T.*(gamma(i)-1))).*(P2(i)./x(i,:)).^(2./gamma(i)) ...
                .*(1-(P2(i)./x(i,:)).^((gamma(i)-1)./gamma(i))));

    P01_residual = Cd*A_inj(i)*P01(i)*sqrt((2*gamma(i)/(R(i)*T*(gamma(i)-1)))*(P2(i)/P01(i))^(2/gamma(i)) ...
                *(1-(P2(i)/P01(i))^((gamma(i)-1)/gamma(i)))) - m_dot(i);
end

figure

plot(x(1,:)/1e5, y(1,:), 'b', 'LineWidth', 2); hold on
plot(x(1,:)/1e5, ym(1,:), 'b--', 'LineWidth', 1.8);

plot(x(2,:)/1e5, y(2,:), 'r', 'LineWidth', 2);
plot(x(2,:)/1e5, ym(2,:), 'r--', 'LineWidth', 1.8);

plot(P01(1)/1e5, m_dot(1), 'bo', 'MarkerFaceColor','b', 'MarkerSize',7);
plot(P01(2)/1e5, m_dot(2), 'ro', 'MarkerFaceColor','r', 'MarkerSize',7);

grid on
xlabel('Upstream Stagnation Pressure, P_0 (bar)', 'FontSize', 12)
ylabel('Mass Flow Rate (MFR) (kg/s)', 'FontSize', 12)
title('Injector Mass Flow Rate vs Upstream Stagnation Pressure', 'FontSize', 14)

legend('Oxygen MFR(P_0)', 'Oxygen Required MFR', 'Methane MFR(P_0)', ...
    'Methane Required MFR', 'Solved P01 (Oxygen)', 'Solved P01 (Methane)',...
    'Location','northeast')

fprintf('\n');
fprintf('====================================================\n');
fprintf('                P01 Solver Residual                \n');
fprintf('====================================================\n\n');

fprintf('Residual                        : %18.15f \n\n', P01_residual);



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



