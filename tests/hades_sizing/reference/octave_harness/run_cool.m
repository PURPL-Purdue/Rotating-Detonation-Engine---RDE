function OUT = run_cool()
  OUT = struct();
  try
    addpath('/tmp/hades_octave/split/HADES_size_chambCooltime'); HADES_size_chambCooltime_body;
  catch err
    fprintf(2, 'v1 plotting error (after results): %s\n', err.message);
  end
  OUT.v1_n = numel(t); OUT.v1_tCool95 = tCool95; OUT.v1_Twall_end_C = T_wall_celsius(end);
  OUT.v1_Tair_end_C = T_air_celsius(end); OUT.v1_T95_wall_C = T95_wall_celsius; OUT.v1_t_end = t(end);
  OUT.v1_t_first = t(1:6)'; OUT.v1_Twall_first = T_wall(1:6)';
end
