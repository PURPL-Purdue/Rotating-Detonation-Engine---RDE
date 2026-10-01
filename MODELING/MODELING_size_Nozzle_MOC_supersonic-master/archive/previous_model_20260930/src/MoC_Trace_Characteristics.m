function paths=MoC_Trace_Characteristics(r)
% Continuous visualization of characteristic directions from the computed
% field. These traces do not replace the solver nodes/edges. Bounding-region
% traces use its explicitly prescribed uniform state.
g=r.geometry; t=r.triple; L=g.period; H=g.height;
count=45; ds=min(L/700,H/150); inset=ds*0.05;
paths=struct('xy',{},'family',{},'region',{});
for region=2:4
    if region==2
        net=r.mesh.products; n=net.nodes;
        xy=net.seeds(:,1:2)+inset*[cos(g.zeta) sin(g.zeta)];
        x=linspace(g.foot+inset,L-inset,count)';
        bottom=max(0,(x-g.refillStart)*tan(g.zeta))+inset;
        top=g.waveHeight+x*tan(t.slipAngle)-inset;
        xy=[xy;x bottom;x top];
        [points,ia]=unique(n(:,1:2),'rows');
        Ft=scatteredInterpolant(points(:,1),points(:,2),n(ia,3),'linear','nearest');
        Fm=scatteredInterpolant(points(:,1),points(:,2),n(ia,4),'linear','nearest');
    elseif region==3
        x=linspace(inset,min(L,(H-g.waveHeight)/tan(t.slipAngle))-inset,count)';
        xy=[x g.waveHeight+x*tan(t.slipAngle)+inset; ...
            x g.waveHeight+x*tan(t.shockAngle)-inset];
    else
        yy=linspace(g.waveHeight+inset,H-inset,count)';
        xx=linspace(inset,L-inset,count)';
        xy=[repmat(inset,count,1) yy;xx repmat(H-inset,count,1)];
    end
    xy=xy(inside(xy,region),:);
    for family=[1 -1]
        for direction=[1 -1]
            p=xy; active=true(size(p,1),1);
            curves=nan(size(p,1),2,1601); curves(:,:,1)=p;
            for step=1:1600
                ix=find(active); if isempty(ix), break; end
                a=angle(p(ix,:),family,region);
                mid=p(ix,:)+direction*ds/2*[cos(a) sin(a)];
                a=angle(mid,family,region);
                pn=p(ix,:)+direction*ds*[cos(a) sin(a)];
                good=inside(pn,region);
                active(ix(~good))=false;
                p(ix(good),:)=pn(good,:);
                curves(ix(good),:,step+1)=pn(good,:);
            end
            for j=1:size(p,1)
                line=squeeze(curves(j,:,:))'; line=line(all(isfinite(line),2),:);
                if size(line,1)>2
                    paths(end+1)=struct('xy',line,'family',family,'region',region); %#ok<AGROW>
                end
            end
        end
    end
end
    function a=angle(p,family,region)
        if region==2, theta=Ft(p(:,1),p(:,2)); M=Fm(p(:,1),p(:,2));
        elseif region==3, theta=t.slipAngle; M=t.M4;
        else, theta=g.zeta; M=t.M5; end
        a=theta+family*asin(1./M);
        if isscalar(a), a=repmat(a,size(p,1),1); end
    end
    function yes=inside(p,region)
        x=p(:,1); y=p(:,2);
        yes=x>=0 & x<=L & y>=0 & y<=H;
        refill=(x<=g.foot & y<g.waveHeight-x/tan(g.zeta)) | ...
            (x>=g.refillStart & y<=max(0,(x-g.refillStart)*tan(g.zeta)));
        rid=2*ones(size(x)); rid(y>g.waveHeight+x*tan(t.slipAngle))=3;
        rid(y>g.waveHeight+x*tan(t.shockAngle))=4; rid(refill)=1;
        yes=yes & rid==region;
    end
end

