function results = HADES_size_ceaRocket(varargin)
% ceaRocket - Modular MATLAB wrapper for NASA CEA rocket equilibrium analysis
%
% Usage examples:
% r = ceaRocket('ox','O2','fuel','H2','of',2.5,'Pc',1000,'PcUnits','psia','Pe',14.7);
% r = ceaRocket('ox','O2','fuel','H2','phi',1.2,'Pc',1000,'areaRatio',8);
%
% Options for mixture specification (choose one):
% 'of' - Oxidizer/Fuel weight ratio
% 'pctFuel' - Fuel percentage by weight
% 'phi' - Equivalence ratio based on fuel/oxidizer
% 'r' - Valence-based equivalence ratio
%
% Pressure units supported: 'psia', 'atm', 'bar', 'mmHg'
% 
% OS Requirements:
% Windows Device
%
% Optional:
% 'Pe' - exit pressure
% 'areaRatio' - Use instead of exit pressure
% 
% The function returns a struct containing the following fields from the
% CEA output file:
% results.AeAt
% results.cstar
% results.cf
% results.isp
% results.ivac

%% PARSE INPUTS
p = inputParser;
addParameter(p,'ox','',@ischar);
addParameter(p,'fuel','',@ischar);
addParameter(p,'Pc',[],@isnumeric);
addParameter(p,'PcUnits','psia',@ischar);
addParameter(p,'Pe',[],@isnumeric);
addParameter(p,'PeUnits','psia',@ischar);
addParameter(p,'areaRatio',[],@isnumeric);
addParameter(p,'of',[],@isnumeric);
addParameter(p,'pctFuel',[],@isnumeric);
addParameter(p,'phi',[],@isnumeric);
addParameter(p,'r',[],@isnumeric);
defaultCeaExe = fullfile(fileparts(mfilename('fullpath')),'CEA','FCEA2.exe');
addParameter(p,'ceaExe',defaultCeaExe,@ischar);

parse(p,varargin{:});
opts = p.Results;

%% VALIDATE INPUTS
if isempty(opts.ox) || isempty(opts.fuel)
    error('Both ''ox'' and ''fuel'' must be specified.');
end
if isempty(opts.Pc)
    error('Chamber pressure ''Pc'' must be specified.');
end
% Default areaRatio if neither Pe nor areaRatio provided
if isempty(opts.Pe) && isempty(opts.areaRatio)
    opts.areaRatio = 10;
end
if ~isempty(opts.Pe) && ~isempty(opts.areaRatio)
    error('Specify either ''Pe'' or ''areaRatio'', not both.');
end

% Determine mixture type
mixCount = sum([~isempty(opts.of), ~isempty(opts.pctFuel), ~isempty(opts.phi), ~isempty(opts.r)]);
if mixCount == 0
    error('Must specify one mixture type: of, pctFuel, phi, or r.');
elseif mixCount > 1
    error('Specify only one mixture type.');
end
if ~isempty(opts.of), mixType='of'; val=opts.of;
elseif ~isempty(opts.pctFuel), mixType='pctFuel'; val=opts.pctFuel;
elseif ~isempty(opts.phi), mixType='phi'; val=opts.phi;
else, mixType='r'; val=opts.r; end

%% PRESSURE CONVERSIONS
Pc_psia = convertPressure(opts.Pc, opts.PcUnits, 'psia');
if ~isempty(opts.Pe)
    Pe_psia = convertPressure(opts.Pe, opts.PeUnits, 'psia');
else
    Pe_psia = [];
end

%% FILE PATHS
ceaDir = fileparts(opts.ceaExe);
inputName  = 'cea_run';
inputFile  = fullfile(ceaDir,[inputName,'.inp']);
outputFile = fullfile(ceaDir,[inputName,'.out']);

%% WRITE CEA INPUT
fid = fopen(inputFile,'w');
fprintf(fid,'problem\n');
fprintf(fid,'    rocket equilibrium\n');
fprintf(fid,'    p,psia=%f\n',Pc_psia);

% Exit condition
if ~isempty(Pe_psia)
    fprintf(fid,'    pi/p=%f\n',Pc_psia/Pe_psia);
else
    fprintf(fid,'    supar=%f\n',opts.areaRatio);
end

% Mixture specification
switch mixType
    case 'of'
        fprintf(fid,'    o/f=%f\n',val);
    case 'pctFuel'
        fprintf(fid,'    o/f=%f\n',(100-val)/val);
    case 'phi'
        fprintf(fid,'    phi=%f\n',val);
    case 'r'
        fprintf(fid,'    r,eq.ratio=%f\n',val);
end

fprintf(fid,'reac\n');
fprintf(fid,'    oxid %s wt=100\n',opts.ox);
fprintf(fid,'    fuel %s wt=100\n',opts.fuel);

fprintf(fid,'output short\n');
fprintf(fid,'only\n');
fprintf(fid,'end\n');
fclose(fid);

%% RUN CEA
orig = pwd;
cd(ceaDir);
runFile = 'run.txt';
fid = fopen(runFile,'w'); fprintf(fid,'%s\n',inputName); fclose(fid);
cmd = sprintf('"%s" < %s > nul', opts.ceaExe, runFile);
system(cmd);
cd(orig);

%% PARSE OUTPUT
results = parseCEAOutput(outputFile);

% Return pressures in original units
results.Pc = opts.Pc;
results.Pe = opts.Pe;
results.PcUnits = opts.PcUnits;
results.PeUnits = opts.PeUnits;

end

%% PARSER
function res = parseCEAOutput(filename)
txt = fileread(filename);
res = struct();

% Initialize with NaN
res.AeAt  = NaN;
res.cstar = NaN;
res.cf    = NaN;
res.isp   = NaN;
res.ivac  = NaN;

lines = splitlines(txt);
perfStart = find(contains(lines,'PERFORMANCE PARAMETERS'),1);

if isempty(perfStart)
    warning('PERFORMANCE PARAMETERS section not found.');
    return
end

for i = perfStart+1 : perfStart+10
    line = strtrim(lines{i});
    if startsWith(line,'Ae/At')
        nums = regexp(line,'[-+]?\d*\.?\d+[eE]?[-+]?\d*','match');
        if ~isempty(nums), res.AeAt  = str2double(nums{1}); end
    elseif contains(line,'CSTAR')
        nums = regexp(line,'[-+]?\d*\.?\d+[eE]?[-+]?\d*','match');
        if ~isempty(nums), res.cstar = str2double(nums{1}); end
    elseif contains(line,'CF')
        nums = regexp(line,'[-+]?\d*\.?\d+[eE]?[-+]?\d*','match');
        if ~isempty(nums), res.cf    = str2double(nums{1}); end
    elseif contains(line,'Ivac')
        nums = regexp(line,'[-+]?\d*\.?\d+[eE]?[-+]?\d*','match');
        if ~isempty(nums), res.ivac  = str2double(nums{1});
        end
    elseif contains(line,'Isp')
        nums = regexp(line,'[-+]?\d*\.?\d+[eE]?[-+]?\d*','match');
        if ~isempty(nums), res.isp   = str2double(nums{1});
        end
    end
end
end

%% HELPERS
function p_out = convertPressure(p_in,unit_in,unit_out)
unit_in = lower(unit_in); unit_out = lower(unit_out);
switch unit_in
    case 'psia', p_pa = p_in*6894.76;
    case 'atm',  p_pa = p_in*101325;
    case 'bar',  p_pa = p_in*1e5;
    case 'mmhg', p_pa = p_in*133.322;
    otherwise, error('Unsupported pressure unit');
end
switch unit_out
    case 'psia', p_out = p_pa/6894.76;
    case 'atm',  p_out = p_pa/101325;
    case 'bar',  p_out = p_pa/1e5;
    case 'mmhg', p_out = p_pa/133.322;
    otherwise, error('Unsupported pressure unit');
end
end
