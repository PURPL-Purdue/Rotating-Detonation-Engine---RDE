function P01 = HADES_size_P_plenum_fanno(Pc, m_dot, A_inj, L, D, T, Cd, gamma, R, mu_ref, T_ref, S_ref, epsilon, M_G)
%HADES_SIZE_P_PLENUM_FANNO  Find plenum stagnation pressure for H2 injector
%
%   Works backwards: Pc (chamber) → Fanno flow (orifice) → Cd contraction → P01 (plenum)
%
%   INPUTS:
%       Pc        - Chamber pressure (Pa)
%       m_dot     - H2 mass flow rate (kg/s)
%       A_inj     - Injector orifice area (m^2)
%       L         - Orifice tube length (mm) 
%       D         - Orifice tube diameter (mm)
%       T         - Static inlet temperature (K)
%       Cd        - Discharge coefficient (-)
%       gamma     - Specific heat ratio for H2 (-)
%       R         - Gas constant for H2 (J/kg-K)
%       mu_ref    - Reference viscosity for H2 (Pa-s)
%       T_ref     - Reference temperature for Sutherland (K)
%       S_ref     - Sutherland constant for H2 (K)
%       epsilon   - Surface roughness (mm)
%       M_G       - Initial Mach number guess (subsonic, e.g. 0.5)
%
%   OUTPUT:
%       P01       - Required plenum stagnation pressure (Pa)

%% ---------- Iteration Setup ----------
tol      = 1e-6;
max_iter = 100;
err      = 1;
iter     = 0;

P2     = Pc * 1.5;   % Initial downstream pressure guess (Pa)
T_stag = T;          % Initial stagnation temperature guess (K)

while err > tol && iter < max_iter

    iter = iter + 1;

    % ---- Sutherland viscosity ----
    mu = mu_ref * (T_ref + S_ref) / (T_stag + S_ref) * (T / T_ref)^(3/2);

    % ---- Reynolds number ----
    Re = (P2 * M_G * sqrt(gamma * R * T) * (D/1000)) / (R * T * mu);

    % ---- Colebrook friction factor ----
    colebrook = @(ff) 1/sqrt(ff) + 2*log10(epsilon/(3.7*D) + 2.51/(Re*sqrt(ff)));
    f = fsolve(colebrook, 0.02, optimoptions('fsolve', 'Display', 'off'));

    % ---- Fanno: solve for inlet Mach number ----
    %   Fanno parameter: fL/D = (1-M^2)/(γM^2) + (γ+1)/(2γ) * ln((γ+1)M^2 / (2+(γ-1)M^2))
    %   We know fL/D from geometry, solve for M_inlet
    fanno_fun = @(M) (1 - M^2)/(gamma*M^2) + ((gamma+1)/(2*gamma)) ...
        * log(((gamma+1)*M^2) / (2 + (gamma-1)*M^2)) - f*L/D;

    M_inlet = fsolve(fanno_fun, M_G, optimoptions('fsolve', 'Display', 'off'));
    M_G = M_inlet;
    if M_inlet <= 0
        error('HADES: Fanno solver returned nonphysical Mach number.');
    end

    % ---- Update stagnation temperature ----
    T_stag_new = T * (1 + ((gamma-1)/2) * M_inlet^2);

    % ---- Fanno: get static pressure at orifice inlet from Pc ----
    %   Isentropic relation from Mach at inlet to choked condition
    P2_new = Pc * (1/M_inlet) * (1/sqrt((2/(gamma+1)) * (1 + ((gamma-1)/2)*M_inlet^2)));

    % ---- Convergence check ----
    err = abs((P2_new - P2) / P2);

    % ---- Update ----
    P2     = P2_new;
    T_stag = T_stag_new;

end

if iter == max_iter
    warning('HADES: Fanno iteration did not fully converge.');
end

%% ---------- Cd contraction — solve for P01 (plenum) ----------
%   Subcritical compressible orifice equation rearranged for P01:
%   m_dot = Cd * A * P01 * sqrt( (2γ / (R*T*(γ-1))) * (P2/P01)^(2/γ) * (1 - (P2/P01)^((γ-1)/γ)) )

P01_fun = @(P0) Cd * A_inj * P0 * sqrt( (2*gamma / (R*T*(gamma-1))) ...
    * (P2/P0)^(2/gamma) * (1 - (P2/P0)^((gamma-1)/gamma)) ) - m_dot;

a = P2 * 1.01;
b = a + 1e8; % 1000 bar search window

if P01_fun(a) * P01_fun(b) < 0
    P01 = fzero(P01_fun, [a b]);
else
    error('HADES: P01 solver could not find a root in search window. Check inputs.');
end

%% Output
fprintf('\n');
fprintf('HADES P_PLENUM FANNO - H2 RESULTS\n');
fprintf("Darcy friction factor: %.4f\n", f);
fprintf('Orifice Inlet Mach Number : %.4f\n',      M_inlet);
fprintf('Orifice Inlet Static Pressure : %.4f bar\n',  P2/1e5);
fprintf('Required Plenum Pressure (P01) : %.4f bar\n',  P01/1e5);
fprintf('Required Plenum Pressure (P01) : %.4f psi\n',  P01*14.504/1e5);
fprintf('Pressure Drop (P01 - Pc) : %.4f bar\n',  (P01-Pc)/1e5);

end
