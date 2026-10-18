function run_original(script, srcdir, outfile, subs_from, subs_to, varnames, postcmd)
% RUN_ORIGINAL  Run an ORIGINAL MATLAB script of the MoC project under GNU Octave
% and save selected workspace variables, for the Python regression tests.
%
%   script    : path of the original .m script (e.g. MOC_2D_steady_irrotational_main.m)
%   srcdir    : folder with the original solver functions (src2/ or src/)
%   outfile   : .mat file written with save -v7
%   subs_from : cell array of regular expressions applied to the script text
%   subs_to   : cell array of replacements (input-section overrides only)
%   varnames  : cell array of workspace variables to save
%   postcmd   : (optional) Octave code evaluated after the script, before saving
%
% Only the *input section* of the script is changed (regexprep); the leading
% clear/clc/close all lines are removed so this function keeps its own state,
% addpath('./src/') is redirected to srcdir, figures are invisible.
  txt = fileread(script);
  txt = regexprep(txt, '^\s*clear\s*;?\s*$', '', 'lineanchors');
  txt = regexprep(txt, '^\s*clc\s*;?\s*$', '', 'lineanchors');
  txt = regexprep(txt, '^\s*close all\s*;?\s*$', '', 'lineanchors');
  txt = regexprep(txt, 'addpath\([^)]*\);', '');
  for k = 1:numel(subs_from)
    new = regexprep(txt, subs_from{k}, subs_to{k}, 'lineanchors');
    if strcmp(new, txt)
      error('run_original: substitution %d (%s) did not match', k, subs_from{k});
    end
    txt = new;
  end
  tmpdir = tempname(); mkdir(tmpdir);
  tmpscript = fullfile(tmpdir, 'moc_case_script.m');
  fid = fopen(tmpscript, 'w'); fwrite(fid, txt); fclose(fid);
  addpath(srcdir);
  set(0, 'defaultfigurevisible', 'off');
  run_original__outfile = outfile; %#ok<NASGU>
  run_original__varnames = varnames; %#ok<NASGU>
  if nargin < 7, postcmd = ''; end
  run_original__postcmd = postcmd; %#ok<NASGU>
  run(tmpscript);
  eval(run_original__postcmd);
  S = struct();
  for k = 1:numel(run_original__varnames)
    S.(run_original__varnames{k}) = eval(run_original__varnames{k});
  end
  save('-v7', run_original__outfile, '-struct', 'S');
  close all;
  rmpath(srcdir);
end
