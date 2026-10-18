function OUT = run_cool3(which)
  OUT = struct();
  try
    if which == 2, addpath('/tmp/hades_octave/split/HADES_size_chambCooltime_v2'); HADES_size_chambCooltime_v2_body; rmpath('/tmp/hades_octave/split/HADES_size_chambCooltime_v2'); else, addpath('/tmp/hades_octave/split/HADES_size_chambCooltime_v3'); HADES_size_chambCooltime_v3_body; end
  catch err
    fprintf(2, 'v%d plotting error (after results): %s\n', which, err.message);
  end
  OUT.t95 = [t_95_base, t_95_N2, t_95_rough, t_95_combo];
  OUT.t90 = [t_90_base, t_90_N2, t_90_rough, t_90_combo];
  OUT.n = [numel(t_base), numel(t_N2), numel(t_rough), numel(t_combo)];
  OUT.Twall_end = [Y_base(end,2), Y_N2(end,2), Y_rough(end,2), Y_combo(end,2)];
end
