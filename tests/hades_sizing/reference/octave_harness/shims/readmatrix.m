function M = readmatrix(fname, varargin)
  % Octave shim: reads the CSV export of an .xlsx (numeric rows only, header skipped)
  [d, n, e] = fileparts(fname);
  csv = fullfile(d, [n '.csv']);
  M = dlmread(csv, ',');
end
