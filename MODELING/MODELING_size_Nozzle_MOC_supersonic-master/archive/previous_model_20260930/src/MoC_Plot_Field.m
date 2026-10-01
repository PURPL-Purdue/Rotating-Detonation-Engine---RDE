function [fig,netfig]=MoC_Plot_Field(r)
% Local gas Mach and detonation propagation Mach are distinct quantities.
f=r.field; g=r.geometry; L=g.period;
paths=MoC_Trace_Characteristics(r);
fig=figure('Name','Unwrapped rotating detonation flow','Color','w', ...
    'Position',[60 60 1350 850]);
layout=tiledlayout(fig,2,2,'TileSpacing','compact');
values={f.Mlab,f.P/1e6,f.T};
titles={'Gas Mach: laboratory frame','Static pressure [MPa, logarithmic color]', ...
    'Static temperature [K]','Continuous wave-fixed characteristic net'};
% Phase is a display translation of one wave pitch, not a periodic solve.
xx=mod(f.x(1,:)+g.displayShift,L); [xx,order]=sort(xx);
for k=1:4
    ax=nexttile(layout); hold(ax,'on');
    if k<4
        surface(ax,repmat(xx,size(f.x,1),1)*1e3,f.y(:,order)*1e3, ...
            zeros(size(f.x)),values{k}(:,order),'EdgeColor','none','FaceColor','flat');
        view(ax,2); colormap(ax,turbo); colorbar(ax);
        if k==2, ax.ColorScale='log'; end
    else
        drawNet(ax);
    end
    boundaries(ax); title(ax,titles{k}); finish(ax);
end
title(layout,sprintf('Straight shock/slip | D/a_{reactants} = %.2f | downstream CJ gas M_{wave} = %.2f', ...
    g.detonationMach,r.inputs.numerics.seedMach));
netfig=figure('Name','Full characteristic net','Color','w','Position',[80 80 1250 650]);
ax=axes(netfig); hold(ax,'on'); drawNet(ax); boundaries(ax); finish(ax);
title(ax,sprintf('C+ blue, C- red | detonation tilt %.2f deg from axial | D/a_{reactants} = %.2f', ...
    rad2deg(g.zeta),g.detonationMach));
    function drawNet(ax)
        for j=1:numel(paths)
            if paths(j).family==1, color=[0.15 0.25 0.85]; else, color=[0.85 0.15 0.2]; end
            line=wrap(paths(j).xy);
            plot(ax,line(:,1)*1e3,line(:,2)*1e3,'Color',color,'LineWidth',0.45);
        end
        x=mod((g.refillStart+L)/2+g.displayShift,L);
        text(ax,x*1e3,g.waveHeight*0.22e3,sprintf('Refill M_{lab} = %.3f',r.injection.M), ...
            'HorizontalAlignment','center','FontSize',8,'BackgroundColor','w');
    end
    function boundaries(ax)
        draw([g.foot 0;0 g.waveHeight],[0.8 0 0],2);
        x=linspace(0,L,500)';
        draw([x g.waveHeight+x*tan(r.triple.slipAngle)],[0 0.5 0],1.4);
        draw([x g.waveHeight+x*tan(r.triple.shockAngle)],[0 0.5 0],1.4);
        x=linspace(g.refillStart,L,100)';
        draw([x (x-g.refillStart)*tan(g.zeta)],[0 0.5 0],1.4);
        function draw(p,color,width)
            % Densify the short angled detonation before phase wrapping.
            if size(p,1)==2, p=p(1,:)+linspace(0,1,40)'.*(p(2,:)-p(1,:)); end
            p=wrap(p); plot(ax,p(:,1)*1e3,p(:,2)*1e3,'Color',color,'LineWidth',width);
        end
    end
    function p=wrap(p)
        p(:,1)=mod(p(:,1)+g.displayShift,L);
        jump=find(abs(diff(p(:,1)))>L/2);
        for j=numel(jump):-1:1
            p=[p(1:jump(j),:);nan(1,2);p(jump(j)+1:end,:)];
        end
    end
    function finish(ax)
        axis(ax,[0 L*1e3 0 g.height*1e3]);
        xlabel(ax,'x - tangential [mm]'); ylabel(ax,'y - axial [mm]');
        ax.Layer='top'; box(ax,'on');
    end
end
