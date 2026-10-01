function f=MoC_Interpolate_Field(mesh,inj,geom,t,opt)
% Region-wise interpolation, with no interpolation through discontinuities.
% region: 1 refill; 2 products; 3 shocked products; 4 bounding products.
[f.x,f.y]=meshgrid(linspace(0,geom.period,opt.gridSize(1)), ...
    linspace(0,geom.height,opt.gridSize(2)));
x=f.x; y=f.y; z=geom.zeta;
local=[x(:) y(:)]*geom.marchRotation-geom.marchOrigin*geom.marchRotation;
sx=reshape(local(:,1),size(x)); sy=reshape(local(:,2),size(x));
slip=interp1(geom.curves(:,1),geom.curves(:,2),sx,'linear','extrap');
shock=interp1(geom.curves(:,1),geom.curves(:,3),sx,'linear','extrap');
refill=max(0,(x-geom.refillStart)*tan(z));
det=geom.waveHeight-x/tan(z);
f.region=2*ones(size(x));
f.region(sy>slip)=3; f.region(sy>shock)=4;
f.region((x<=geom.foot & y<det) | (x>=geom.refillStart & y<=refill))=1;
f.theta=nan(size(x)); f.Mwave=f.theta; f.T=f.theta; f.P=f.theta;
f.covered=false(size(x)); f.analyticClosure=false(size(x));
% Product interpolation only within its characteristic hull.
n=mesh.products.nodes;
[nxy,ia]=unique(n(:,1:2),'rows');
F=scatteredInterpolant(nxy(:,1),nxy(:,2),n(ia,3),'linear','none');
idx=f.region==2; f.theta(idx)=F(x(idx),y(idx));
F.Values=n(ia,4); f.Mwave(idx)=F(x(idx),y(idx));
f.covered(idx)=isfinite(f.Mwave(idx));
f=thermo(f,idx,t.gas);
% Independent post-shock mesh, including streamline-transported shock loss.
idx=f.region==3; n=mesh.shocked.nodes;
[nxy,ia]=unique(n(:,1:2),'rows');
F=scatteredInterpolant(nxy(:,1),nxy(:,2),n(ia,3),'linear','none');
f.theta(idx)=F(x(idx),y(idx)); F.Values=n(ia,4); f.Mwave(idx)=F(x(idx),y(idx));
F.Values=log(n(ia,5)); p0=exp(F(x(idx),y(idx)));
f.T(idx)=t.gas.T0./(1+(t.gas.gamma-1)*f.Mwave(idx).^2/2);
f.P(idx)=p0.*(f.T(idx)/t.gas.T0).^(t.gas.gamma/(t.gas.gamma-1));
% The constant post-shock state is an exact extension under the restored
% straight shock / uniform upstream closure, not arbitrary extrapolation.
missing=idx & ~isfinite(f.P);
if mesh.shocked.uniformClosure
    f.theta(missing)=t.slipAngle; f.Mwave(missing)=t.M4;
    f=thermo(f,missing,t.shockGas);
    f.analyticClosure(missing)=true;
end
f.covered(idx)=isfinite(f.P(idx));
idx=f.region==4;
f.theta(idx)=z; f.Mwave(idx)=t.M5;
f=thermo(f,idx,t.gas); f.covered(idx)=true; f.analyticClosure(idx)=true;
R=t.gas.R*ones(size(x)); gamma=t.gas.gamma*ones(size(x));
idx=f.region==1;
f.P(idx)=inj.P; f.T(idx)=inj.T; R(idx)=inj.R; gamma(idx)=inj.gamma;
f.theta(idx)=z; f.Mwave(idx)=hypot(geom.waveSpeedX,inj.V)/sqrt(inj.gamma*inj.R*inj.T);
f.covered(idx)=true; f.analyticClosure(idx)=true;
f.rho=f.P./(R.*f.T); f.a=sqrt(gamma.*R.*f.T);
f.uWave=f.Mwave.*f.a.*cos(f.theta); f.vWave=f.Mwave.*f.a.*sin(f.theta);
f.uLab=f.uWave-geom.waveSpeedX; f.vLab=f.vWave;
f.Mlab=hypot(f.uLab,f.vLab)./f.a;
f.nu=MoC_PM(f.Mwave,t.gas.gamma);
% Refill is updated directly, without supersonic compatibility coordinates,
% even though the Galilean wave-frame Mach there is generally supersonic.
f.nu(idx)=NaN;
f.Kplus=f.theta-f.nu; f.Kminus=f.theta+f.nu;
f.gamma=gamma; f.R=R;
f.unresolved=~f.covered;
if any(f.unresolved,'all')
    warning('MoC:Coverage','%.2f%% of grid is outside the characteristic hull (left NaN).', ...
        100*mean(f.unresolved,'all'));
end
end
function f=thermo(f,idx,gas)
f.T(idx)=gas.T0./(1+(gas.gamma-1)*f.Mwave(idx).^2/2);
f.P(idx)=gas.P0*(f.T(idx)/gas.T0).^(gas.gamma/(gas.gamma-1));
end
