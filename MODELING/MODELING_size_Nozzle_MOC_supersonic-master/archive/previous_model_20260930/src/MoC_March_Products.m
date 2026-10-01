function mesh=MoC_March_Products(geom,t,opt)
% Semi-Lagrangian MOC: trace C+/C- feet onto each preceding initial-value
% line, interpolate invariants, then apply midpoint direction correction.
% This remeshes the old src2 predictor/corrector net without centered fans.
% March along the wave normal to avoid the near-vertical CJ characteristics.
z=geom.zeta; rot=[cos(z) -sin(z);sin(z) cos(z)]; foot=[geom.foot 0];
g=t.gas.gamma; n=opt.seedCount; h=geom.waveHeight/cos(z);
knee=([geom.refillStart 0]-foot)*rot;
xend=([geom.period geom.height]-foot)*rot(:,1);
x=0; y=linspace(0,h,n)'; theta=zeros(n,1); M=repmat(t.seedMach,n,1);
nodes=[zeros(n,1) y theta M]; edges=zeros(0,3); slipPressureMismatch=0;
for step=1:opt.maxRows
    nu=MoC_PM(M,g); kp=theta-nu; km=theta+nu;
    slopes=[tan(theta+asin(1./M)) tan(theta-asin(1./M))];
    if any(~isfinite(slopes),'all') || any(abs(theta)+asin(1./M)>=pi/2)
        error('MoC:March','A characteristic reverses in the marching frame.');
    end
    dx=min(0.65*min(diff(y))/max(abs(slopes),[],'all'),xend-x);
    if x<knee(1) && x+dx>knee(1), dx=knee(1)-x; end
    if dx<=eps(max(1,x)), error('MoC:March','Characteristic step collapsed.'); end
    xn=x+dx;
    if xn<=knee(1), low=-xn*tan(z); lowAngle=-z;
    else, low=knee(2); lowAngle=0; end
    high=h+xn*tan(t.delta);
    yn=linspace(low,high,n)';
    th=interp1(y,theta,yn,'linear','extrap');
    mach=interp1(y,M,yn,'linear','extrap');
    mach=max(mach,1+1e-8);
    for iteration=1:20
        sp=0.5*(interp1(y,slopes(:,1),min(max(yn,y(1)),y(end)))+tan(th+asin(1./mach)));
        sm=0.5*(interp1(y,slopes(:,2),min(max(yn,y(1)),y(end)))+tan(th-asin(1./mach)));
        yp=yn-dx*sp; ym=yn-dx*sm;
        plus=interp1(y,kp,min(max(yp,y(1)),y(end)),'linear');
        minus=interp1(y,km,min(max(ym,y(1)),y(end)),'linear');
        % Incoming family + tangency determines the outgoing boundary family.
        plus(yp<y(1))=2*lowAngle-minus(1);
        minus(ym>y(end))=2*t.delta-plus(end);
        plus(1)=2*lowAngle-minus(1);
        minus(end)=2*t.delta-plus(end);
        newTheta=(plus+minus)/2; newNu=(minus-plus)/2;
        if any(newNu<=0), error('MoC:Sonic','Product march reached a sonic boundary.'); end
        newM=MoC_PM(newNu,g,true);
        change=max(abs(newTheta-th)+abs(newM-mach));
        th=newTheta; mach=newM;
        if change<1e-8, break; end
    end
    if change>1e-5, error('MoC:Corrector','Characteristic corrector did not converge.'); end
    new=[repmat(xn,n,1) yn th mach];
    upperGlobal=new(end,1:2)*rot'+foot;
    if upperGlobal(2)<=geom.height
        pressure=t.gas.P0/(1+(g-1)*mach(end)^2/2)^(g/(g-1));
        slipPressureMismatch=max(slipPressureMismatch,abs(pressure-t.P)/t.P);
    end
    newIds=size(nodes,1)+(1:n)'; nodes=[nodes;new]; %#ok<AGROW>
    % Explicit interpolated characteristic feet, rather than connecting to
    % the nearest mesh node (which would misrepresent characteristic slopes).
    for family=[1 -1]
        if family==1, yf=yp; else, yf=ym; end
        xf=repmat(x,n,1);
        crossedLow=yf<y(1); crossedHigh=yf>y(end);
        for j=find(crossedLow | crossedHigh)'
            if crossedLow(j), b0=y(1); b1=low; else, b0=y(end); b1=high; end
            fraction=(yf(j)-b0)/(yf(j)-yn(j)+b1-b0);
            fraction=min(max(fraction,0),1);
            xf(j)=x+dx*fraction;
            yf(j)=b0+(b1-b0)*fraction;
        end
        footTheta=interp1(y,theta,min(max(yf,y(1)),y(end)));
        footM=interp1(y,M,min(max(yf,y(1)),y(end)));
        for j=find(crossedLow | crossedHigh)'
            fraction=(xf(j)-x)/dx;
            if crossedLow(j), k=1; else, k=n; end
            footTheta(j)=(1-fraction)*theta(k)+fraction*th(k);
            footM(j)=MoC_PM((1-fraction)*MoC_PM(M(k),g)+ ...
                fraction*MoC_PM(mach(k),g),g,true);
        end
        footIds=size(nodes,1)+(1:n)';
        nodes=[nodes;xf yf footTheta footM]; %#ok<AGROW>
        edges=[edges;footIds newIds repmat(family,n,1)]; %#ok<AGROW>
    end
    x=xn; y=yn; theta=th; M=mach;
    if x>=xend-1e-12, break; end
end
if x<xend-1e-10
    error('MoC:Extent','maxRows exhausted before covering chamber; increase maxRows.');
end
nodes(:,1:2)=nodes(:,1:2)*rot'+foot; nodes(:,3)=nodes(:,3)+z;
mesh=struct('nodes',nodes,'edges',edges,'gas',t.gas,'rows',step, ...
    'seeds',nodes(1:n,:),'slipPressureMismatch',slipPressureMismatch);
end
