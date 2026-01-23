function results = HADES_size_ceaDet(varargin)
% ceaDet - MATLAB wrapper for NASA CEA detonation equilibrium analysis
%
% Usage examples:
% d = HADES_size_ceaDet('ox','O2','fuel','H2','of',2.5,'P0',1000,'P0Units','psia','T0',300,'T0Units','K');
%
% Options for mixture specification (choose one):
% 'of' - Oxidizer/Fuel weight ratio
% 'pctFuel' - Fuel percentage by weight
% 'phi' - Equivalence ratio based on fuel/oxidizer
% 'r' - Valence-based equivalence ratio
%
% Pressure units: psia, atm, bar, mmHg
% Temperature units: K, C, F, R
%
% OS Requirement:
% Windows Device
%
% Output structure:
% results.cjVel - Chapman-Jouguet velocity
% results.detMach - Detonation Mach Number
% results.P_ratio - P/P1
% results.T_ratio - T/T1
% results.M_ratio - M/M1
% results.RHO_ratio - RHO/RHO1
% results.P_burned_bar - Burned Pressure
% results.T_cj - Burned Gas Temperature


%% PARSE INPUTS
p = inputParser;
addParameter(p,'ox','',@ischar);
addParameter(p,'fuel','',@ischar);
addParameter(p,'P0',[],@isnumeric);
addParameter(p,'P0Units','psia',@ischar);
addParameter(p,'T0',[],@isnumeric);
addParameter(p,'T0Units','K',@ischar);
addParameter(p,'of',[],@isnumeric);
addParameter(p,'pctFuel',[],@isnumeric);
addParameter(p,'phi',[],@isnumeric);
addParameter(p,'r',[],@isnumeric);
addParameter(p,'ceaExe', fullfile(pwd,'CEA','FCEA2.exe'), @ischar);
parse(p,varargin{:});
opts = p.Results;

%% VALIDATE INPUTS
if isempty(opts.ox) || isempty(opts.fuel)
    error('Both oxidizer and fuel must be specified.');
end
if isempty(opts.P0)
    error('Initial pressure P0 must be specified.');
end
if isempty(opts.T0)
    error('Initial temperature T0 must be specified.');
end

mixCount = sum([~isempty(opts.of), ~isempty(opts.pctFuel), ~isempty(opts.phi), ~isempty(opts.r)]);
if mixCount ~= 1
    error('Specify exactly one mixture type.');
end

if ~isempty(opts.of), mixType='of'; val=opts.of;
elseif ~isempty(opts.pctFuel), mixType='pctFuel'; val=opts.pctFuel;
elseif ~isempty(opts.phi), mixType='phi'; val=opts.phi;
else, mixType='r'; val=opts.r;
end

%% UNIT CONVERSIONS
P0_psia = convertPressure(opts.P0, opts.P0Units, 'psia');
T0_K    = convertTemperature(opts.T0, opts.T0Units, 'K');

%% FILE PATHS
ceaDir = fileparts(opts.ceaExe);
inputName = 'cea_det';
inputFile = fullfile(ceaDir,[inputName,'.inp']);
outputFile = fullfile(ceaDir,[inputName,'.out']);

%% WRITE CEA INPUT
fid = fopen(inputFile,'w');
fprintf(fid,'problem\n');
fprintf(fid,'    det\n');
fprintf(fid,'    p,psia=%f\n',P0_psia);
fprintf(fid,'    t,k=%f\n',T0_K);

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
results = parseCEAOutputDet(outputFile);

results.P0 = opts.P0;
results.T0 = opts.T0;
results.P0Units = opts.P0Units;
results.T0Units = opts.T0Units;

end

%% DETONATION PARSER
function res = parseCEAOutputDet(filename)

txt = fileread(filename);
lines = splitlines(txt);

res = struct( ...
    'P_ratio',NaN, ...
    'T_ratio',NaN, ...
    'M_ratio',NaN, ...
    'rho_ratio',NaN, ...
    'detMach',NaN, ...
    'cjVel',NaN, ...
    'P_burned_bar',NaN, ...
    'T_cj',NaN);

inBurnedGas = false;

for i=1:length(lines)
    L = strtrim(lines{i});

    % detonation parameters
    if startsWith(L,'P/P1')
        res.P_ratio = sscanf(L,'P/P1 %f');
    elseif startsWith(L,'T/T1')
        res.T_ratio = sscanf(L,'T/T1 %f');
    elseif startsWith(L,'M/M1')
        res.M_ratio = sscanf(L,'M/M1 %f');
    elseif startsWith(L,'RHO/RHO1')
        res.rho_ratio = sscanf(L,'RHO/RHO1 %f');
    elseif contains(L,'DET MACH NUMBER')
        nums = regexp(L,'[-+]?\d*\.?\d+','match');
        res.detMach = str2double(nums{1});
    elseif contains(L,'DET VEL')
        nums = regexp(L,'[-+]?\d*\.?\d+','match');
        res.cjVel = str2double(nums{1});
    end

    % Burned gas section flag
    if contains(L,'BURNED GAS')
        inBurnedGas = true;
    end

    if inBurnedGas
        if startsWith(L,'P, BAR')
            nums = regexp(L,'[-+]?\d*\.?\d+','match');
            res.P_burned_bar = str2double(nums{1});
        elseif startsWith(L,'T, K')
            nums = regexp(L,'[-+]?\d*\.?\d+','match');
            res.T_cj = str2double(nums{1});
        end
    end
end
end

%% PRESSURE CONVERSION
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

%% TEMPERATURE CONVERSION
function T_out = convertTemperature(T_in,unit_in,unit_out)

unit_in = upper(unit_in); unit_out = upper(unit_out);

switch unit_in
    case 'K', T_K = T_in;
    case 'C', T_K = T_in + 273.15;
    case 'F', T_K = (T_in-32)*5/9 + 273.15;
    case 'R', T_K = T_in*5/9;
    otherwise, error('Unsupported temperature unit');
end

switch unit_out
    case 'K', T_out = T_K;
    case 'C', T_out = T_K - 273.15;
    case 'F', T_out = (T_K-273.15)*9/5 + 32;
    case 'R', T_out = T_K*9/5;
end
end
