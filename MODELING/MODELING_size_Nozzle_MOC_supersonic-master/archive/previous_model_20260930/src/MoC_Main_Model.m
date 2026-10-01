function result = MoC_Main_Model(options)
%MOC_MAIN_MODEL Straight-boundary, wave-fixed unwrapped RDC model.
% r = MoC_Main_Model(struct('numerics',struct('plot',false)));
% Solver units: SI. Input geometry/port dimensions: mm. See MODEL_NOTES.md.
if nargin == 0, options = struct; end
c.chemistry = struct('ox','O2','fuel','CH4','mixture_spec','phi', ...
    'mixture_spec_value',1.3,'useCEA',true);
c.chamber = struct('outer_chamber_diameter',52.243, ...
    'inner_chamber_diameter',37.167,'chamber_length',29.21);
% P3/T3 are STATIC port conditions. Cp is computed from gamma and R.
c.injection = struct('injector_area',3.7138,'N_holes',36, ...
    'P3',1e6,'T3',398,'mdot',0.59,'Pa',1e6,'R',315,'gamma1',1.3569);
% Optional round-port diameter overrides injector_area (one port per set).
c.injection.port_diameter = [];
% Manual chemistry is an explicit alternative, never a silent CEA fallback.
c.products = struct('cjVel',2563,'P',8e6,'T',3905,'R',430.47,'gamma',1.1405);
c.numerics = struct('seedCount',31,'maxRows',7000,'seedMach',1.02, ...
    'gridSize',[301 151],'plot',true,'save',true,'waveCount',1, ...
    'refillFraction',0.65,'displayPhase',0.58);
c = mergeOptions(c,options);
modelDir = fileparts(mfilename('fullpath'));
validateattributes(c.numerics.seedCount,{'numeric'},{'scalar','integer','>=',5});
validateattributes(c.numerics.maxRows,{'numeric'},{'scalar','integer','positive'});
validateattributes(c.numerics.seedMach,{'numeric'},{'scalar','>',1});
validateattributes(c.numerics.waveCount,{'numeric'},{'scalar','integer','positive'});
validateattributes(c.numerics.refillFraction,{'numeric'},{'scalar','>',0,'<',1});
validateattributes(c.numerics.gridSize,{'numeric'},{'vector','numel',2,'integer','>=',3});
[circumference,height] = MoC_Calculate_Domain_Size(c.chamber);
inj = MoC_Injection_Velocity(c.injection,c.chamber);
if c.chemistry.useCEA
    chem = c.chemistry;
    cea = HADES_size_ceaDet('ox',chem.ox,'fuel',chem.fuel, ...
        chem.mixture_spec,chem.mixture_spec_value,'P0',inj.P/1e5, ...
        'P0Units','bar','T0',inj.T,'T0Units','K','ceaExe',getCEAPath());
    products = struct('cjVel',cea.cjVel,'P',cea.P_burned_bar*1e5, ...
        'T',cea.T_cj,'R',cea.R_specific,'gamma',cea.gamma_burned);
else
    cea = struct; products = c.products;
end
vals = [products.cjVel products.P products.T products.R products.gamma-1];
if any(~isfinite(vals) | vals<=0)
    error('MoC:Chemistry','Product properties must be finite and physical.');
end
if inj.V >= products.cjVel
    error('MoC:Triangle','Reactant speed must be below CJ wave speed.');
end
geometry.period = circumference/c.numerics.waveCount;
geometry.height = height;
geometry.zeta = asin(inj.V/products.cjVel); % Thesis Eq. 2.21
geometry.waveSpeedX = products.cjVel*cos(geometry.zeta);
geometry.detonationMach = products.cjVel/sqrt(inj.gamma*inj.R*inj.T);
geometry.displayShift = c.numerics.displayPhase*geometry.period;
validateattributes(c.numerics.displayPhase,{'numeric'},{'scalar','finite','>=',0,'<',1});
% Explicit refill-duration closure, not a converged injector-blocking model.
geometry.waveHeight = c.numerics.refillFraction*geometry.period* ...
    sin(geometry.zeta)*cos(geometry.zeta);
if geometry.waveHeight >= height
    error('MoC:Geometry','Wave height exceeds chamber length; reduce refillFraction.');
end
geometry.foot = geometry.waveHeight*tan(geometry.zeta);
geometry.refillStart = geometry.period-geometry.waveHeight/tan(geometry.zeta);
triple = MoC_Resolve_Triple_Point(products,inj,geometry,c.numerics.seedMach);
[mesh,geometry] = MoC_Solve_Field(products,inj,geometry,triple,c.numerics);
field = MoC_Interpolate_Field(mesh,inj,geometry,triple,c.numerics);
result = struct('inputs',c,'injection',inj,'chemistry',cea, ...
    'products',products,'geometry',geometry,'triple',triple,'mesh',mesh,'field',field);
result.diagnostics = struct('productNodes',size(mesh.products.nodes,1), ...
    'shockNodes',size(mesh.shocked.nodes,1), ...
    'productCoverage',mean(field.covered(field.region==2)), ...
    'pressureMatchRelative',triple.pressureResidual, ...
    'model','Straight shock/slip; triple-point match only; prescribed bounding gas');
result.diagnostics.slipPressureMismatch = mesh.products.slipPressureMismatch;
result.diagnostics.seedOffset = c.numerics.seedMach-1;
if c.numerics.plot
    [result.figure,result.netFigure] = MoC_Plot_Field(result);
end
if c.numerics.save
    outputDir = fullfile(modelDir,'..','results');
    if ~isfolder(outputDir), mkdir(outputDir); end
    saved = result;
    if isfield(saved,'figure'), saved = rmfield(saved,{'figure','netFigure'}); end
    save(fullfile(outputDir,'MoC_result.mat'),'saved');
    if c.numerics.plot
        exportgraphics(result.figure,fullfile(outputDir,'unwrapped_flow.png'),'Resolution',180);
        exportgraphics(result.netFigure,fullfile(outputDir,'characteristic_net.png'),'Resolution',180);
    end
end
fprintf('Injection %.1f m/s (M_lab %.3f); zeta %.2f deg; wave height %.2f mm\n', ...
    inj.V,inj.M,rad2deg(geometry.zeta),1e3*geometry.waveHeight);
fprintf('MOC nodes: products %d, shocked %d; product grid coverage %.1f%%\n', ...
    result.diagnostics.productNodes,result.diagnostics.shockNodes, ...
    100*result.diagnostics.productCoverage);
fprintf('Wave propagation D/a_reactants = %.3f; downstream CJ seed M_wave = %.3f\n', ...
    geometry.detonationMach,c.numerics.seedMach);
end
function a = mergeOptions(a,b)
if ~isstruct(b) || ~isscalar(b), error('MoC:Options','Options must be a scalar structure.'); end
names = fieldnames(b);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(a,name), error('MoC:Options','Unknown option: %s',name); end
    if isstruct(a.(name)), a.(name) = mergeOptions(a.(name),b.(name));
    else, a.(name) = b.(name); end
end
end
