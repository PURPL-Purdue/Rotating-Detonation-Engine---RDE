function [mesh,geom] = MoC_Solve_Field(~,~,geom,t,opt)
% Line-seeded characteristic mesh. The predictor/corrector construction and
% characteristic intersections follow src2/..._internal_point.m. We integrate
% its planar compatibility equation as theta-/+nu instead of solving for u,v.
% Nodes: [x y theta Mach], edges: [parent child family], angles radians.
mesh.products=MoC_March_Products(geom,t,opt);
% Separate characteristic net propagated inward from interpolated shock seeds.
% Under the uniform bounding-state closure the straight shock produces a
% uniform state. It is an exact constant-state test for the second MOC net.
n=opt.seedCount;
xend=min(geom.period,(geom.height-geom.waveHeight)/tan(t.shockAngle));
s=linspace(0,1,n)';
shockEnds=[0 geom.waveHeight; xend geom.waveHeight+xend*tan(t.shockAngle)];
xy=interp1([0;1],shockEnds,s);
row=[xy repmat([t.slipAngle t.M4],n,1)];
snodes=row; sedges=zeros(0,3); ids=(1:n)';
for generation=1:n-1
    next=zeros(size(row,1)-1,4); nextIds=zeros(size(next,1),1);
    for j=1:size(next,1)
        % C+ from the upper seed and C- from the lower seed: inward normal.
        next(j,:)=MoC_Interior_Point(row(j+1,:),row(j,:),t.gas.gamma);
        snodes(end+1,:)=next(j,:); nextIds(j)=size(snodes,1); %#ok<AGROW>
        sedges(end+1,:)=[ids(j+1) nextIds(j) 1]; %#ok<AGROW>
        sedges(end+1,:)=[ids(j) nextIds(j) -1]; %#ok<AGROW>
    end
    row=next; ids=nextIds;
end
% Keep only points inside the shock/slip wedge. Add exact characteristic
% intersections at the slip, preserving its uniform tangential state.
for j=2:n
    b=boundary(snodes(j,:),-1,[0 geom.waveHeight],t.slipAngle,t.gas.gamma);
    snodes(end+1,:)=b; sedges(end+1,:)=[j size(snodes,1) -1]; %#ok<AGROW>
end
inside=snodes(:,1)>=-1e-10 & snodes(:,1)<=geom.period & ...
    snodes(:,2)>=geom.waveHeight+snodes(:,1)*tan(t.slipAngle)-1e-10 & ...
    snodes(:,2)<=geom.waveHeight+snodes(:,1)*tan(t.shockAngle)+1e-10;
keep=inside(sedges(:,1)) & inside(sedges(:,2));
map=cumsum(inside); sedges=sedges(keep,:);
sedges(:,1:2)=map(sedges(:,1:2));
mesh.shocked=struct('nodes',snodes(inside,:),'edges',sedges,'gas',t.shockGas, ...
    'seeds',snodes(1:n,:),'uniformClosure',true);
end

function b=boundary(a,family,origin,angle,gamma)
nu=MoC_PM(a(4),gamma);
invariant=a(3)-family*nu;
newNu=family*(angle-invariant);
if newNu<=0, error('MoC:Boundary','Boundary condition reaches sonic/subsonic flow.'); end
M=MoC_PM(newNu,gamma,true);
phi=0.5*(a(3)+family*asin(1/a(4))+angle+family*asin(1/M));
d=[cos(phi);sin(phi)]; wall=[cos(angle);sin(angle)];
A=[d -wall];
if rcond(A)<1e-12, error('MoC:Boundary','Characteristic tangent to boundary.'); end
q=A\(origin-a(1:2))';
xy=a(1:2)+q(1)*d';
b=[xy angle M];
end

