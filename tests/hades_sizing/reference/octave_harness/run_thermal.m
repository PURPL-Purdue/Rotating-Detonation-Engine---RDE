function OUT = run_thermal()
  OUT = struct();
  try
    addpath('/tmp/hades_octave/split/HADES_size_thermal_v1'); HADES_size_thermal_v1_body;
  catch err
    fprintf(2, 'thermal plotting error (after results): %s\n', err.message);
  end
  OUT.time_to_limit = time_to_limit'; OUT.h_inner = h_inner; OUT.T_rec_g = T_rec_g; OUT.Re_g = Re_g;
  OUT.Nu_g = Nu_g; OUT.T_inner = T_inner'; OUT.T_end_profile = T(end,:);
end
