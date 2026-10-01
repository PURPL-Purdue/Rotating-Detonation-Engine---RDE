function exePath = getCEAPath()
% Find FCEA2.exe next to the source files, independent of the current folder.
sourceDir = fileparts(mfilename('fullpath'));
ceaDir = fullfile(sourceDir, 'CEA');
exePath = fullfile(ceaDir, 'FCEA2.exe');

if ~isfile(exePath)
    error('MoC:CEA', 'CEA executable not found: %s', exePath);
end
end
