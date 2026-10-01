function node=MoC_Interior_Point(plus,minus,gamma)
% Planar compatibility integration and midpoint characteristic correction.
% Adapted from src2/MOC_2D_steady_irrotational_internal_point.m; theta=atan2
% convention and absolute convergence/conditioning checks replace old code.
kp=plus(3)-MoC_PM(plus(4),gamma);
km=minus(3)+MoC_PM(minus(4),gamma);
theta=0.5*(kp+km); nu=0.5*(km-kp);
if nu<=0, error('MoC:Sonic','Characteristic intersection is not supersonic.'); end
M=MoC_PM(nu,gamma,true);
% Correct endpoint directions (trapezoidal characteristic approximation).
p=0.5*(plus(3)+asin(1/plus(4))+theta+asin(1/M));
m=0.5*(minus(3)-asin(1/minus(4))+theta-asin(1/M));
dp=[cos(p);sin(p)]; dm=[cos(m);sin(m)]; A=[dp -dm];
if rcond(A)<1e-12, error('MoC:Parallel','Characteristic lines are nearly parallel.'); end
q=A\(minus(1:2)-plus(1:2))';
xy=plus(1:2)+q(1)*dp';
node=[xy theta M];
end
