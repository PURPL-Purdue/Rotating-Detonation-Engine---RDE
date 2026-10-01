function report=run_tests()
% Meaningful numerical checks; no Optimization Toolbox is required.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','src'));
g=1.4;
assert(isnan(MoC_PM(0.5,g)));
assert(abs(MoC_PM(2,g)-0.460413682082695)<1e-12);
assert(max(abs(MoC_PM(MoC_PM([1 1.02 2 5],g),g,true)-[1 1.02 2 5]))<1e-10);
% Uniform flow: exact intersection and compatibility preservation.
a=[0 0 0 2]; b=[0 1 0 2]; q=MoC_Interior_Point(a,b,g);
assert(norm(q-[sqrt(3)/2 0.5 0 2])<1e-10);
a=[0 0 0.03 2]; b=[0 1 0.07 2.2]; q=MoC_Interior_Point(a,b,g);
assert(abs(q(3)-MoC_PM(q(4),g)-(a(3)-MoC_PM(a(4),g)))<1e-10);
assert(abs(q(3)+MoC_PM(q(4),g)-(b(3)+MoC_PM(b(4),g)))<1e-10);
base=struct('chemistry',struct('useCEA',false), ...
    'numerics',struct('seedCount',11,'plot',false,'save',false,'gridSize',[121 61]));
r=MoC_Main_Model(base); f=r.field;
assert(all(f.covered,'all') && all(isfinite(f.P),'all'));
assert(all(f.P>0 & f.T>0 & f.rho>0,'all'));
assert(max(abs(f.P-f.rho.*f.R.*f.T),[],'all')/max(f.P,[],'all')<1e-12);
assert(abs(r.injection.massResidual)<1e-12 && abs(r.injection.energyResidual)<1e-7);
refill=f.region==1;
assert(max(abs(f.uLab(refill)))<1e-8);
assert(max(abs(f.vLab(refill)-r.injection.V))<1e-8);
assert(all(isnan(f.nu(refill))));
assert(all(f.Mlab(refill)<1) && all(f.Mwave(refill)>1));
assert(r.triple.pressureResidual<1e-10);
assert(isfinite(r.diagnostics.slipPressureMismatch) && r.diagnostics.slipPressureMismatch>=0);
assert(abs(r.geometry.detonationMach-r.products.cjVel/sqrt(r.injection.gamma*r.injection.R*r.injection.T))<1e-12);
assert(max(r.mesh.shocked.nodes(:,5))-min(r.mesh.shocked.nodes(:,5))<1e-9);
assert(all(r.mesh.shocked.nodes(:,5)<=r.triple.gas.P0*(1+1e-8)));
assert(r.triple.shockGas.P0<r.triple.gas.P0);
assert(abs(r.triple.shockGas.T0-r.triple.gas.T0)<1e-9);
% Angled detonation and straight shock/slip geometry.
assert(r.geometry.foot>0);
assert(abs(r.geometry.foot/r.geometry.waveHeight-tan(r.geometry.zeta))<1e-12);
s=r.mesh.shocked.nodes;
q=s(:,1:2)*r.geometry.marchRotation-r.geometry.marchOrigin*r.geometry.marchRotation;
lo=interp1(r.geometry.curves(:,1),r.geometry.curves(:,2),q(:,1),'linear','extrap');
hi=interp1(r.geometry.curves(:,1),r.geometry.curves(:,3),q(:,1),'linear','extrap');
assert(all(q(:,2)>=lo-1e-8 & q(:,2)<=hi+1e-8));
% Independent CEA parsing regression against known chemistry, not live files.
c=HADES_size_ceaDet('ox','O2','fuel','CH4','phi',1.3,'P0',2.5, ...
    'P0Units','bar','T0',283,'T0Units','K','ceaExe',getCEAPath());
assert(abs(c.cjVel-2564.7)<1);
assert(abs(c.R_specific-8314.462618/19.286)<0.1);
assert(abs(c.gamma_unburned-1.3568)<0.001 && abs(c.gamma_burned-1.1406)<0.001);
% Invalid geometry and incomplete marches must fail explicitly.
expectError(@() MoC_Main_Model(struct('chamber',struct('inner_chamber_diameter',60))),'MoC:Geometry');
bad=base; bad.numerics.seedMach=1;
try, MoC_Main_Model(bad); error('test:missingError','Sonic seed accepted');
catch e, assert(~strcmp(e.identifier,'test:missingError')); end
bad=base; bad.numerics.maxRows=1;
expectError(@() MoC_Main_Model(bad),'MoC:Extent');
% Refinement is reported, not misrepresented as external validation.
finer=base; finer.numerics.seedCount=21; rr=MoC_Main_Model(finer);
mask=f.region==2 & rr.field.region==2;
report.pressureRefinementL1=mean(abs(f.P(mask)-rr.field.P(mask)))/mean(rr.field.P(mask));
assert(report.pressureRefinementL1<0.15,'Mesh refinement change exceeded 15%%.');
report.slipPressureMismatch=rr.diagnostics.slipPressureMismatch;
report.status='PASS'; disp(report);
end
function expectError(fun,id)
try, fun(); error('test:missingError','Expected failure did not occur');
catch e, assert(strcmp(e.identifier,id),'Unexpected error: %s',e.message); end
end



