set(0, 'defaultfigurevisible', 'off');
addpath('/tmp/hades_octave/shims');
C = struct();
C.v1 = run_cool(); C.v2 = run_cool3(2); C.v3 = run_cool3(3); C.thermal = run_thermal();
fid = fopen('/tmp/hades_octave/octave_cool.json', 'w'); fprintf(fid, '%s', jsonencode(C)); fclose(fid);
