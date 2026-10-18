set(0, 'defaultfigurevisible', 'off');
addpath('/tmp/hades_octave/shims');
OUT = struct();
% --- convert
cv = {1,'atm','bar'; 1,'bar','atm'; 100,'psia','bar'; 2,'bar','mmHg'; 10,'m/s','ft/s'; 3,'ft/s','m/s'; ...
      300,'K','C'; 25,'degC','F'; 500,'R','K'; 77,'F','K'; 0.25,'in','m'; 2,'ft','in'; 1,'lbms','kgs'; 3,'kgs','lbms'};
for i = 1:size(cv,1), OUT.convert(i) = HADES_size_convert(cv{i,1}, cv{i,2}, cv{i,3}); end
% --- throat
OUT.throat = [HADES_size_throat(0.93, 9.8e5, 2330, 1.22, 333), HADES_size_throat(1, 2e6, 300, 1.4, 287)];
% --- pdrop_fanno
[a,b,c,d] = HADES_size_pdrop_fanno(20e5, 4e-3, 0.1, 3.2e-6, 0.01, 1.4, 1.8e-5); OUT.pdrop1 = [a b c d];
[a,b,c,d] = HADES_size_pdrop_fanno(40e5, 6e-3, 0.3, 1e-5, 0.02, 1.4, 8.8e-6); OUT.pdrop2 = [a b c d];
% --- hoop stress, other inputs
[a,b] = HADES_size_HoopStressTemps(2, 2.5, 0.101325, 30.7, 1.5); OUT.hoop1 = [a b];
[a,b] = HADES_size_HoopStressTemps(1.375, 2, 0.101325, 45, 2.0); OUT.hoop2 = [a b];
% --- plenum fanno, other inputs
OUT.plenum2 = HADES_size_P_plenum_fanno(8e5, 5e-4, pi*0.25*0.0015^2, 8, 1.5, 290, 0.65, 1.4, 4124.2, 84e-7, 273, 72, 3.2e-3, 0.4);
% --- annulus gap / fill height, other inputs
[a,b,c] = HADES_size_annulusgap_fillheight(1.1, 0.25, 2, 1980, 1100, 30, 3, 420); OUT.fill2 = [a b c];
% --- chamber pressure, other inputs
OUT.chp2 = HADES_size_chamberPressure(2, 10.0, 1950.3, 15.2, 30);
% --- injectors etc.
OUT.h2 = HADES_size_h2Inj(0.06, 150); OUT.air = zeros(1,2); [OUT.air(1), OUT.air(2)] = HADES_size_airInj(1.8, 8, 2.2, 290);
OUT.geom = HADES_size_geometry(2.5, 0.5, 8); OUT.thrust = HADES_size_thrust(1500, 2.1); OUT.P0 = HADES_size_P0_v1();
OUT.flux = HADES_size_max_flux_range(2.1, 2000);
% --- sonic nozzle sweep (script)
[OUT.sonic_p, OUT.sonic_dia_mm] = run_sonic();
% --- cell_size2 (superseded script): fit-based cell width at phi=0.9
[OUT.cell_size2_cellwidth, OUT.cell_size2_coeffs] = run_cell_size2();
fid = fopen('/tmp/hades_octave/octave_units.json', 'w'); fprintf(fid, '%s', jsonencode(OUT)); fclose(fid);
