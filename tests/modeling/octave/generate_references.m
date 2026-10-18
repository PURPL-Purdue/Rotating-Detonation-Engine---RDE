function generate_references(repo_root, outdir)
% GENERATE_REFERENCES  Build the Octave reference data used by tests/modeling.
%
%   octave --no-gui --norc --eval "addpath([pwd '/tests/modeling/octave']); generate_references(pwd, 'tests/modeling/data')"
%
% Every case runs the ORIGINAL MATLAB sources (MODELING_size_*/...) unmodified,
% except for input-section overrides applied to a temporary copy of the driver
% script (see run_original.m).  Results are saved as .mat (v7).
%
% Source location: the originals were removed from the working tree after the
% port.  Recover them from git history (see MODELING/moc/PORTING_NOTES.md,
% "Regenerating the Octave references"):
%   mkdir -p /tmp/rde_src && git archive 2a30c95 | tar -x -C /tmp/rde_src
%   export RDE_MATLAB_SOURCE_ROOT=/tmp/rde_src
% The folder holding MODELING_size_Nozzle_MOC_supersonic-master and
% MODELING_size_MoC_unchanged_backup is $RDE_MATLAB_SOURCE_ROOT/MODELING (a root
% that already is that folder is accepted too), falling back to <repo_root>/MODELING.
  addpath(fileparts(mfilename('fullpath')));  % absolute: run() changes directory
  if ~exist(outdir, 'dir'), mkdir(outdir); end
  outdir = make_absolute_filename(outdir);
  repo_root = make_absolute_filename(repo_root);
  src_root = getenv('RDE_MATLAB_SOURCE_ROOT');
  if isempty(src_root)
    src_root = repo_root;
  end
  src_root = make_absolute_filename(src_root);
  if exist(fullfile(src_root, 'MODELING'), 'dir') || ~exist(fullfile(src_root, 'MODELING_size_Nozzle_MOC_supersonic-master'), 'dir')
    src_root = fullfile(src_root, 'MODELING');
  end
  team = fullfile(src_root, 'MODELING_size_Nozzle_MOC_supersonic-master');
  upstream = fullfile(src_root, 'MODELING_size_MoC_unchanged_backup');
  if ~exist(team, 'dir') || ~exist(upstream, 'dir')
    error(['generate_references: original MATLAB sources not found under %s. Recover them with\n' ...
           '  mkdir -p /tmp/rde_src && git archive 2a30c95 | tar -x -C /tmp/rde_src\n' ...
           'and set RDE_MATLAB_SOURCE_ROOT=/tmp/rde_src'], src_root);
  end
  printf('original sources: %s\n', src_root);
  rde = fullfile(team, 'RDE_Toolbox_2');

  nozzle_vars = {'X','Y','U','V','LENG_INDI','indI','ind_patch_tri','ind_patch_quad', ...
                 'patch_i_tri','patch_j_tri','patch_i_quad','patch_j_quad','LOC_PATCHES', ...
                 'xvnull','uvnull','ythroat','Mach','pressure','wallnodesJ', ...
                 'chocked_mach','chocked_pressure','case_geom','case_params'};
  post = ['[Mach,pressure] = MOC_2D_steady_irrotational_get_thermo(U,V,params);' ...
          'wallnodesJ = geom.NI+1 : size(X,2); wallnodesJ = wallnodesJ(X(1,wallnodesJ) < geom.xe);' ...
          'chocked = get_Laval_theory(geom,params,X(1,wallnodesJ),Y(1,wallnodesJ));' ...
          'chocked_mach = chocked.mach; chocked_pressure = chocked.pressure;' ...
          'case_geom = geom; case_params = params;'];
  set_field = @(name, value) {sprintf('^(\\s*%s\\s*=)[^;%%]*;', regexprep(name, '\.', '\\.')), ...
                              sprintf('$1 %s ;', value)};

  cases = {
    % name                 script                                   src                         overrides
    'moc_team',            fullfile(team,'MOC_2D_steady_irrotational_main.m'),     fullfile(team,'src2'),     {};
    'moc_upstream',        fullfile(upstream,'MOC_2D_steady_irrotational_main.m'), fullfile(upstream,'src'),  {};
    'moc_upstream_te0',    fullfile(upstream,'MOC_2D_steady_irrotational_main.m'), fullfile(upstream,'src'),  {{'geom.te','0'}};
    'moc_team_te5',        fullfile(team,'MOC_2D_steady_irrotational_main.m'),     fullfile(team,'src2'),     {{'geom.te','5'}};
    'moc_team_weakshock',  fullfile(team,'MOC_2D_steady_irrotational_main.m'),     fullfile(team,'src2'),     {{'params.PRatio','0.8'}};
    'moc_team_strongshock',fullfile(team,'MOC_2D_steady_irrotational_main.m'),     fullfile(team,'src2'),     {{'params.PRatio','0.4'}};
    'moc_upstream_weakshock', fullfile(upstream,'MOC_2D_steady_irrotational_main.m'), fullfile(upstream,'src'), {{'params.PRatio','0.8'}};
    'moc_upstream_NI21',   fullfile(upstream,'MOC_2D_steady_irrotational_main.m'), fullfile(upstream,'src'),  {{'geom.NI','21'}};
  };
  for c = 1:rows(cases)
    ov = cases{c,4};
    from = {}; to = {};
    for k = 1:numel(ov)
      st = set_field(ov{k}{1}, ov{k}{2});
      from{end+1} = st{1}; to{end+1} = st{2};
    end
    printf('--- %s\n', cases{c,1});
    try
      run_original(cases{c,2}, cases{c,3}, fullfile(outdir, [cases{c,1} '.mat']), from, to, nozzle_vars, post);
    catch err
      % the original code stops with an error for this case: record the message
      error_message = err.message; %#ok<NASGU>
      save('-v7', fullfile(outdir, [cases{c,1} '.mat']), 'error_message');
      printf('    original code raised: %s\n', err.message);
      close all;
    end
  end

  % ---- stand-alone scripts of src2 -------------------------------------------
  printf('--- minimum_length_nozzle\n');
  run_original(fullfile(team,'src2','minimum_length_nozzle.m'), fullfile(team,'src2'), ...
               fullfile(outdir,'mln.mat'), {}, {}, ...
               {'abs','ord','theta','nu','mu','Mach','Kminus','Kplus','Wall','absWall','ordWall', ...
                'absp','ordp','Machp','area_ratio','theta_max','nu_out'});
  printf('--- expansion\n');
  run_original(fullfile(team,'src2','expansion.m'), fullfile(team,'src2'), ...
               fullfile(outdir,'expansion.mat'), {}, {}, {'delta','nu1','mu1','M2','p1static','p2static'});
  printf('--- oblique_shock\n');
  run_original(fullfile(team,'src2','oblique_shock.m'), fullfile(team,'src2'), ...
               fullfile(outdir,'oblique_shock.mat'), {}, {}, {'machmax','betadev'});

  % ---- RDE toolbox functions -------------------------------------------------
  printf('--- RDE toolbox\n');
  shim = tempname(); mkdir(shim);
  % Octave has no optimoptions: map the fsolve options used by MoC_Resolve_Triple_Point
  fid = fopen(fullfile(shim,'optimoptions.m'),'w');
  fprintf(fid, ['function o = optimoptions(solver, varargin)\n' ...
                '  map = struct(''Display'',''Display'',''FunctionTolerance'',''TolFun'',''StepTolerance'',''TolX'');\n' ...
                '  args = {};\n  for k = 1:2:numel(varargin)\n    args(end+1:end+2) = {map.(varargin{k}), varargin{k+1}};\n  end\n' ...
                '  o = optimset(args{:});\nend\n']);
  fclose(fid);
  % MoC_Injection_Velocity.m does not parse (fprintf(..., out.);): run a copy with the
  % 5th printed value replaced by 0 so the computed fields can be compared.
  txt = fileread(fullfile(rde,'src','MoC_Injection_Velocity.m'));
  txt = strrep(txt, 'out.T4(k), out.);', 'out.T4(k), 0);');
  fid = fopen(fullfile(shim,'MoC_Injection_Velocity.m'),'w'); fwrite(fid, txt); fclose(fid);
  addpath(shim);
  addpath(fullfile(rde,'src'));
  addpath(shim);  % shim first (MoC_Injection_Velocity copy)

  ChamberDim.outer_chamber_diameter = 52.243;
  ChamberDim.inner_chamber_diameter = 37.167;
  ChamberDim.chamber_length         = 29.21;
  InjV.injector_area = 3.7138;
  InjV.N_holes = 36;
  iTripleParam.P1 = 250000; iTripleParam.P2 = 8e6; iTripleParam.T2 = 3905;
  iTripleParam.Vcj = 2563; iTripleParam.gamma2 = 1.1405; iTripleParam.R1 = 430.47;
  InjV.P3 = 1e6; InjV.T3 = 398; InjV.mdot = 0.59; InjV.gamma1 = 1.3569;
  InjV.Pa = 1e6; InjV.R = 315; InjV.Cp = 315;

  out = MoC_Injection_Velocity(InjV, ChamberDim);
  [beta, delta, P_match, M3, M1_p] = MoC_Resolve_Triple_Point(iTripleParam);
  [xMax, yMax] = MoC_Calculate_Domain_Size(ChamberDim);
  inj_V4 = out.V4; inj_P4 = out.P4; inj_T4 = out.T4; inj_Vwave = out.Vwave;
  save('-v7', fullfile(outdir,'rde.mat'), 'inj_V4','inj_P4','inj_T4','inj_Vwave', ...
       'beta','delta','P_match','M3','M1_p','xMax','yMax');
  rmpath(shim); rmpath(fullfile(rde,'src'));
  printf('done\n');
end
