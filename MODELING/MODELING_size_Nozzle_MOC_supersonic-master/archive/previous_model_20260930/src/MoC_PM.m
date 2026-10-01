function value = MoC_PM(value,gamma,inverse)
% Prandtl-Meyer compatibility coordinate in radians, NOT a velocity potential.
% No fan generation. Undefined for M<1.
if nargin<3, inverse=false; end
if inverse
    limit = pi/2*(sqrt((gamma+1)/(gamma-1))-1);
    if any(~isfinite(value(:)) | value(:)<0 | value(:)>=limit)
        error('MoC:Characteristic','Incompatible characteristic invariants.');
    end
    lo=ones(size(value)); hi=2*lo;
    while any(MoC_PM(hi,gamma)<value,'all'), hi=hi*2; end
    for k=1:48
        mid=(lo+hi)/2; below=MoC_PM(mid,gamma)<value;
        lo(below)=mid(below); hi(~below)=mid(~below);
    end
    value=(lo+hi)/2;
else
    M=value; value=nan(size(M)); valid=isfinite(M) & M>=1;
    b=sqrt(M(valid).^2-1);
    value(valid)=sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1))*b)-atan(b);
end
end
