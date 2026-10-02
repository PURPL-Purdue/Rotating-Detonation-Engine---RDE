function run_tests()
% Check only the live CEA interface. The physical model is not implemented.
sourceDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
addpath(sourceDir);

originalFolder = pwd;
outputFolder = fullfile(sourceDir, '..', 'results', 'CEA_test');

cea = HADES_size_ceaDet( ...
    'ox', 'O2', ...
    'fuel', 'CH4', ...
    'phi', 1.3, ...
    'P0', 2.5, ...
    'P0Units', 'bar', ...
    'T0', 283, ...
    'T0Units', 'K', ...
    'ceaExe', getCEAPath(), ...
    'outputDir', outputFolder);

assert(strcmp(pwd, originalFolder), 'CEA changed the caller working folder.');
assert(isfile(cea.inputFile), 'CEA input was not retained.');
assert(isfile(cea.outputFile), 'CEA output was not retained.');
assert(contains(fileread(cea.inputFile), 'phi=1.300000'));
assert(contains(fileread(cea.outputFile), 'DETONATION PROPERTIES'));
assert(abs(cea.cjVel - 2564.7) < 1, 'Unexpected CEA velocity or parsing error.');
assert(abs(cea.R_specific - 8314.462618 / 19.286) < 0.1);
assert(abs(cea.gamma_unburned - 1.3568) < 0.001);
assert(abs(cea.gamma_burned - 1.1406) < 0.001);

P0SweepStart_bar = 2.5;
P0SweepEnd_bar = 25.5;
nSweepPoints = 33;
P0Sweep_bar = linspace(P0SweepStart_bar, P0SweepEnd_bar, nSweepPoints);
cjVelSweep = nan(size(P0Sweep_bar));

for k = 1:numel(P0Sweep_bar)
    ceaSweep = HADES_size_ceaDet( ...
        'ox', 'O2', ...
        'fuel', 'CH4', ...
        'phi', 1.3, ...
        'P0', P0Sweep_bar(k), ...
        'P0Units', 'bar', ...
        'T0', 283, ...
        'T0Units', 'K', ...
        'ceaExe', getCEAPath());
    cjVelSweep(k) = ceaSweep.cjVel;
end

assert(all(isfinite(cjVelSweep)), 'CEA returned an invalid CJ velocity in the pressure sweep.');
figure('Name', 'CJ velocity vs. initial pressure');
plot(P0Sweep_bar, cjVelSweep, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Initial pressure P_0 (bar)');
ylabel('CJ velocity (m/s)');
title('CJ velocity vs. initial pressure');

fprintf('CEA interface check passed. No RDE model calculation was tested.\n');
end
