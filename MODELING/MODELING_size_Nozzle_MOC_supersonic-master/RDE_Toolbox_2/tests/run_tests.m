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

fprintf('CEA interface check passed. No RDE model calculation was tested.\n');
end
