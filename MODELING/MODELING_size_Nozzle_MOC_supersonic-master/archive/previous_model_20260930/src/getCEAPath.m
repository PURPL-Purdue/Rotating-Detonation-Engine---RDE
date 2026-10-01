function path = getCEAPath()
% Locate bundled CEA deterministically; never open an interactive file dialog.
path=fullfile(fileparts(mfilename('fullpath')),'CEA','FCEA2.exe');
if ~isfile(path)
    error('MoC:CEA','Bundled CEA executable missing: %s. Set chemistry.useCEA=false for manual products.',path);
end
end
