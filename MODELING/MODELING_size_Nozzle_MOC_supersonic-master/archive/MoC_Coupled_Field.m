function [mesh,geom]=MoC_Coupled_Field(geom,t,opt)
% Rotational, steady 2-D MOC with a fitted contact and weak oblique shock.
% Characteristic compatibility: dtheta +/- A dp = 0,
% A=sqrt(M^2-1)/(gamma*p*M^2). Total enthalpy is constant; total pressure
% is transported on streamlines, retaining the spatially varying shock loss.
% Midpoint compatibility and semi-Lagrangian interpolation of C+/C- feet
% extend the planar velocity compatibility used by the archived src2 solver.
g=t.gas.gamma; z=geom.zeta; rot=[cos(z) -sin(z);sin(z) cos(z)];
origin=[geom.foot 0]; h=geom.waveHeight/cos(z); n=opt.seedCount;
knee=([geom.refillStart 0]-origin)*rot;
xend=([geom.period geom.height]-origin)*rot(:,1);
% A finite first wedge avoids coincident slip/shock grid lines at the triple
% point. Its length decreases with refinement; the detonation seeds remain
% on the actual wave line at x=0.
x0=h/(n-1)*0.05; x=x0;
slip=h+x*tan(t.delta); shock=h+x*tan(t.beta);
y1=linspace(-x*tan(z),slip,n)'; y2=linspace(slip,shock,n)';
pseed=t.gas.P0/(1+(g-1)*t.seedMach^2/2)^(g/(g-1));
s1=[zeros(n,1) repmat(pseed,n,1) repmat(t.gas.P0,n,1)];
s1(1,1)=-z; s1(end,:)=[t.delta t.P t.gas.P0];
s2=repmat([t.delta t.P t.shockGas.P0],n,1);
seedY=linspace(0,h,n)';
nodes1=[zeros(n,1) seedY zeros(n,1) repmat(t.seedMach,n,1) repmat(t.gas.P0,n,1)];
seedX=linspace(0,x0,n)';
nodes2=[seedX h+seedX*tan(t.beta) repmat(t.delta,n,1) ...
    repmat(t.M4,n,1) repmat(t.shockGas.P0,n,1)];
edges1=zeros(0,3); edges2=zeros(0,3);
curve=[0 h h t.delta t.beta t.P t.shockGas.P0];
slipAngle=t.delta; beta=t.beta; worstPressure=0; worstAngle=0;
for step=1:opt.maxRows
    [m1,~]=state(s1); [m2,~]=state(s2);
    slope1=[tan(s1(:,1)+asin(1./m1)) tan(s1(:,1)-asin(1./m1))];
    slope2=[tan(s2(:,1)+asin(1./m2)) tan(s2(:,1)-asin(1./m2))];
    if any(abs(s1(:,1))+asin(1./m1)>=pi/2) || any(abs(s2(:,1))+asin(1./m2)>=pi/2)
        error('MoC:March','A characteristic turns upstream in the marching frame.');
    end
    dx=min([0.45*min(diff(y1))/max(abs(slope1),[],'all'), ...
        0.45*min(diff(y2))/max(abs(slope2),[],'all'),xend-x]);
    if x<knee(1) && x+dx>knee(1), dx=knee(1)-x; end
    if dx<1e-14, error('MoC:March','Coupled step collapsed.'); end
    xn=x+dx;
    if xn<=knee(1), low=-xn*tan(z); lowAngle=-z;
    else, low=knee(2); lowAngle=0; end
    slipNew=slip+dx*tan(slipAngle); shockNew=shock+dx*tan(beta);
    yy1=linspace(low,slipNew,n)'; yy2=linspace(slipNew,shockNew,n)';
    a=interp1(y1,s1,min(max(yy1,y1(1)),y1(end)));
    b=interp1(y2,s2,min(max(yy2,y2(1)),y2(end)));
    for iteration=1:35
        [fp1,fm1]=feet(y1,s1,yy1,a,dx);
        [fp2,fm2]=feet(y2,s2,yy2,b,dx);
        % Entropy/total pressure advection along the streamline family.
        yf=yy2-dx*tan(0.5*(b(:,1)+interp1(y2,s2(:,1),min(max(yy2,y2(1)),y2(end)))));
        p02=interp1(y2,s2(:,3),min(max(yf,y2(1)),y2(end)));
        crossed=yf>y2(end);
        p02(crossed)=b(end,3);
        p02(1)=t.shockGas.P0; % entropy of the contact streamline
        aa=interior(fp1,fm1,repmat(t.gas.P0,n,1));
        bb=interior(fp2,fm2,p02);
        % Common pressure and direction from the two incoming families.
        contact=interior(fp1(end,:),fm2(1,:),[t.gas.P0 t.shockGas.P0]);
        aa(end,:)=[contact(1) contact(2) t.gas.P0];
        bb(1,:)=[contact(1) contact(2) t.shockGas.P0];
        aa(1,:)=wall(fm1(1,:),lowAngle,t.gas.P0);
        [bb(end,:),newBeta]=shockBoundary(fp2(end,:));
        err=max([abs(aa(:,1)-a(:,1));abs(bb(:,1)-b(:,1)); ...
            abs(aa(:,2)-a(:,2))/t.gas.P0;abs(bb(:,2)-b(:,2))/t.gas.P0]);
        a=aa; b=bb;
        if err<2e-7, break; end
    end
    if err>2e-5
        error('MoC:Corrector','Coupled characteristic corrector failed at x=%g (error %g).',xn,err);
    end
    worstPressure=max(worstPressure,abs(a(end,2)-b(1,2))/a(end,2));
    worstAngle=max(worstAngle,abs(a(end,1)-b(1,1)));
    [nodes1,edges1]=record(nodes1,edges1,xn,yy1,a,x,y1,s1,dx);
    [nodes2,edges2]=record(nodes2,edges2,xn,yy2,b,x,y2,s2,dx);
    x=xn; y1=yy1; y2=yy2; s1=a; s2=b;
    slip=slipNew; shock=shockNew; slipAngle=a(end,1); beta=newBeta;
    curve(end+1,:)=[x slip shock slipAngle beta a(end,2) b(end,3)]; %#ok<AGROW>
    if x>=xend-1e-12, break; end
end
if x<xend-1e-10, error('MoC:Extent','Increase maxRows to cover the coupled field.'); end
nodes1(:,1:2)=nodes1(:,1:2)*rot'+origin; nodes1(:,3)=nodes1(:,3)+z;
nodes2(:,1:2)=nodes2(:,1:2)*rot'+origin; nodes2(:,3)=nodes2(:,3)+z;
mesh.products=struct('nodes',nodes1,'edges',edges1,'gas',t.gas,'rows',step, ...
    'seeds',nodes1(1:n,:),'slipPressureMismatch',worstPressure);
mesh.shocked=struct('nodes',nodes2,'edges',edges2,'gas',t.shockGas, ...
    'uniformClosure',false,'seeds',nodes2(1:n,:));
geom.curves=curve;
geom.slipXY=curve(:,[1 2])*rot'+origin;
geom.shockXY=curve(:,[1 3])*rot'+origin;
geom.startupLength=x0;
geom.contactAngleResidual=worstAngle;
geom.marchRotation=rot; geom.marchOrigin=origin;
    function [M,A]=state(s)
        ratio=s(:,3)./s(:,2);
        M2=2/(g-1)*(ratio.^((g-1)/g)-1);
        if any(~isfinite(M2) | M2<=1)
            error('MoC:Sonic','Coupled MOC encountered a non-supersonic state.');
        end
        M=sqrt(M2); A=sqrt(M2-1)./(g*s(:,2).*M2);
    end
    function [pfoot,mfoot]=feet(y,s,yn,sn,dx)
        [mn,~]=state(sn); [mo,~]=state(s);
        for fam=[1 -1]
            oldSlope=tan(s(:,1)+fam*asin(1./mo));
            initial=interp1(y,oldSlope,min(max(yn,y(1)),y(end)));
            slope=0.5*(initial+tan(sn(:,1)+fam*asin(1./mn)));
            yf=yn-dx*slope;
            sf=interp1(y,s,min(max(yf,y(1)),y(end)));
            for side=[1 numel(y)]
                if side==1, crossed=yf<y(1); else, crossed=yf>y(end); end
                fraction=(yf(crossed)-y(side))./(yf(crossed)-yn(crossed)+yn(side)-y(side));
                fraction=min(max(fraction,0),1);
                sf(crossed,:)=(1-fraction).*s(side,:)+fraction.*sn(side,:);
            end
            if fam==1, pfoot=sf; else, mfoot=sf; end
        end
    end
    function s=interior(pf,mf,p0)
        % Integrate both compatibility equations with trapezoidal coefficients.
        [~,ap]=state(pf); [~,am]=state(mf);
        if size(p0,2)==2, p0p=p0(1); p0m=p0(2);
        else, p0p=p0; p0m=p0; end
        lo=min(pf(:,2),mf(:,2))*0.02;
        hi=min(p0p,p0m)/(1+(g-1)/2)^(g/(g-1))*(1-1e-9);
        for it=1:42
            p=(lo+hi)/2;
            [~,anp]=state([zeros(size(p)) p p0p]);
            [~,anm]=state([zeros(size(p)) p p0m]);
            thp=pf(:,1)-0.5*(ap+anp).*(p-pf(:,2));
            thm=mf(:,1)+0.5*(am+anm).*(p-mf(:,2));
            above=thp>thm; lo(above)=p(above); hi(~above)=p(~above);
        end
        p=(lo+hi)/2;
        if any(abs(thp-thm)>2e-6)
            error('MoC:Compatibility','No supersonic compatibility root (residual %g).',max(abs(thp-thm)));
        end
        s=[0.5*(thp+thm) p p0p];
    end
    function s=wall(mf,angle,p0)
        [~,am]=state(mf); lo=mf(2)*0.01;
        hi=p0/(1+(g-1)/2)^(g/(g-1))*(1-1e-9);
        for it=1:42
            p=(lo+hi)/2; [~,an]=state([angle p p0]);
            th=mf(1)+0.5*(am+an)*(p-mf(2));
            if th<angle, lo=p; else, hi=p; end
        end
        s=[angle (lo+hi)/2 p0];
        if abs(th-angle)>2e-6
            error('MoC:Compatibility','No supersonic wall compatibility root.');
        end
    end
    function [s,beta]=shockBoundary(pf)
        [~,ap]=state(pf); mu=asin(1/t.M5);
        lo=mu+1e-9; hi=t.betaMax;
        flo=residual(lo); fhi=residual(hi);
        if flo*fhi>0
            error('MoC:Shock','No attached supersonic shock at x=%g; residual [%g %g].',xn,flo,fhi);
        end
        for it=1:35
            mid=(lo+hi)/2; fm=residual(mid);
            if fm*flo>0, lo=mid; flo=fm; else, hi=mid; end
        end
        beta=(lo+hi)/2; s=shockState(beta);
        function f=residual(beta)
            ss=shockState(beta); [~,an]=state(ss);
            f=ss(1)-pf(1)+0.5*(ap+an)*(ss(2)-pf(2));
        end
    end
    function s=shockState(beta)
        mn=t.M5^2*sin(beta)^2;
        delta=atan(2*cot(beta)*(mn-1)/(t.M5^2*(g+cos(2*beta))+2));
        p=t.boundingPressure*(1+2*g/(g+1)*(mn-1));
        m2=((1+(g-1)*mn/2)/(g*mn-(g-1)/2))/sin(beta-delta)^2;
        p0=p*(1+(g-1)*m2/2)^(g/(g-1));
        s=[delta p p0];
    end
    function [nodes,edges]=record(nodes,edges,xn,yn,sn,x,y,s,dx)
        [mn,~]=state(sn); ids=size(nodes,1)+(1:n)';
        nodes=[nodes;repmat(xn,n,1) yn sn(:,1) mn sn(:,3)];
        [pf,mf]=feet(y,s,yn,sn,dx);
        for fam=[1 -1]
            if fam==1, sf=pf; else, sf=mf; end
            [mf0,~]=state(sf);
            slope=0.5*(tan(sf(:,1)+fam*asin(1./mf0))+tan(sn(:,1)+fam*asin(1./mn)));
            yf=yn-dx*slope; xf=repmat(x,n,1);
            for side=[1 n]
                if side==1, crossed=yf<y(1); else, crossed=yf>y(end); end
                frac=(yf(crossed)-y(side))./(yf(crossed)-yn(crossed)+yn(side)-y(side));
                frac=min(max(frac,0),1);
                xf(crossed)=x+dx*frac; yf(crossed)=y(side)+(yn(side)-y(side))*frac;
            end
            fi=size(nodes,1)+(1:n)';
            nodes=[nodes;xf yf sf(:,1) mf0 sf(:,3)]; %#ok<AGROW>
            edges=[edges;fi ids repmat(fam,n,1)]; %#ok<AGROW>
        end
    end
end
