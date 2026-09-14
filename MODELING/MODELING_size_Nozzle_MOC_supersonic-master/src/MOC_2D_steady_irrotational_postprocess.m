function [fig1,fig2,fig3] = MOC_2D_steady_irrotational_postprocess(geom,params,plots,indI,LENG_INDI,...
                                                                   X,Y,U,V,...
                                                                   ind_tri,ind_quad,...
                                                                   i_tri,j_tri,i_quad,j_quad,LOC_PATCHES)
%
% Post-process the results
%
  disp('Postprocessing the results...')
  
  [Mach,pressure,density,temperature,sound] = MOC_2D_steady_irrotational_get_thermo(U,V,params);

  fig1=figure(1);
  plot(X(1,1:indI)/geom.yt,Y(1,1:indI)/geom.yt,'k','linewidth',2); grid on; hold on;
  xlabel('Axial location x / y_t'); ylabel('Radial location y / y_t')

  if (plots.patches>0)
     plot_colours = true;
     switch plots.patches_data
       case 0
          plot_colours = false;
          output_patch = Mach;
       case 1
          output_patch = Mach;
       case 2
          output_patch = pressure;
       case 3
          output_patch = temperature;
       case 4
          output_patch = density;
       case 5
          output_patch = sqrt(U.^2 + V.^2);
       case 6
          output_patch = atand(V./U);
       otherwise
          error('Unknown type of data to display')
     end
     
     indP = 0; % to update patch_x_tri, patch_y_tri, patch_c_tri
     indP2 = 0; % to locate in LOC_PATCHES
     indN1 = 1; % to locate the nodes
     indN2 = 1; % to locate the nodes
     indF1 = 1; % to locate the faces
     indF2 = 1; % to locate the faces
     colours = [0 0 0;1 0 0;0 0 1;1 0 1];

     % Below is the fix for the patches not rendering the colors. Creating
     % a Patch error. Basically, what was happening before is that you had
     % to adjust the plots.patches.xlim so that the number of patches
     % didn't exceed the number of x nodes. Now, they are automatically
     % cutoff and adjusted with the two modifiers below. 
     nodesC_tri  = nan(3*ind_tri ,1);
     nodesC_quad = nan(4*ind_quad,1);
     nodesRGB_tri  = zeros(3*ind_tri ,3);
     nodesRGB_quad = zeros(4*ind_quad,3);

     for I=1:ind_tri
       nodes_tri(indN1:indN1+2,1:2) = [ X(j_tri(1,I),i_tri(1,I))  ,  Y(j_tri(1,I),i_tri(1,I)) ;...
                                        X(j_tri(2,I),i_tri(2,I))  ,  Y(j_tri(2,I),i_tri(2,I)) ;...
                                        X(j_tri(3,I),i_tri(3,I))  ,  Y(j_tri(3,I),i_tri(3,I))   ] ; % [x y]
       if ( I == LOC_PATCHES(indP2+1,1)+1 )
         indP2 = indP2 + 1; % change of region [ initial-value line ; initial expansion ; nozzle ; plume ]
       end
       %% NEW LINE(64) TO FIX CHARACTERISTIC LINES NOT SHOWING
       nodesRGB_tri(indN1:indN1+2,1:3) = repmat(colours(indP2+1,:),3,1);
       if ( sum( nodes_tri(indN1:indN1+2,1)<plots.patches_xlim ) == 3 )
         % Plot only patches with abscissae below the maximum plot_patches_xlim chosen by user
         faces_tri(indF1,1:3) = indN1:indN1+2 ;
         facesC_tri(indF1,1:3) = colours(indP2+1,:);
         nodesC_tri(indN1:indN1+2,1) = [ output_patch(j_tri(1,I),i_tri(1,I)) ; ...
                                         output_patch(j_tri(2,I),i_tri(2,I)) ; ...
                                         output_patch(j_tri(3,I),i_tri(3,I))   ] ;
         indF1 = indF1 + 1;
       end
       indN1 = indN1 + 3;
     end
     
     indP = 0;
     indP2 = 0;
     for I=1:ind_quad
       nodes_quad(indN2:indN2+3,1:2) = [ X(j_quad(1,I),i_quad(1,I))  ,  Y(j_quad(1,I),i_quad(1,I)) ;...
                                         X(j_quad(2,I),i_quad(2,I))  ,  Y(j_quad(2,I),i_quad(2,I)) ;...
                                         X(j_quad(3,I),i_quad(3,I))  ,  Y(j_quad(3,I),i_quad(3,I)) ;...
                                         X(j_quad(4,I),i_quad(4,I))  ,  Y(j_quad(4,I),i_quad(4,I))   ] ; % [x y]
       if ( I == LOC_PATCHES(indP2+1,2)+1 )
         indP2 = indP2 + 1; % change of region [ initial-value line ; initial expansion ; nozzle ; plume ]
       end
       %% NEW LINE(88) FROM CLAUDE TO FIX CHARACTERISTIC LINES NOT SHOWING
       nodesRGB_quad(indN2:indN2+3,1:3) = repmat(colours(indP2+1,:),4,1);
       if ( sum( nodes_quad(indN2:indN2+3,1)<plots.patches_xlim ) == 4 )
         % Plot only patches with abscissae below the maximum plot_patches_xlim chosen by user
         faces_quad(indF2,1:4) =  indN2:indN2+3 ;
         facesC_quad(indF2,1:3) = colours(indP2+1,:);
         nodesC_quad(indN2:indN2+3,1) = [ output_patch(j_quad(1,I),i_quad(1,I)) ; ...
                                          output_patch(j_quad(2,I),i_quad(2,I)) ; ...
                                          output_patch(j_quad(3,I),i_quad(3,I)) ; ...
                                          output_patch(j_quad(4,I),i_quad(4,I))   ] ;
         indF2 = indF2 + 1;
       end
       indN2 = indN2 + 4;
     end
     
     % %% OLD PLOT COLOURS
     % if (plot_colours);
     % % Plot the left- and right-running characteristics with data
     % patch('Faces',faces_tri,'Vertices',nodes_tri,'FaceVertexCData',nodesC_tri,...
     %         'FaceColor','interp','LineWidth',1,'LineStyle','-');
     % patch('Faces',faces_quad,'Vertices',nodes_quad,'FaceVertexCData',nodesC_quad,...
     %         'FaceColor','interp','LineWidth',1,'LineStyle','-');
     % else
     % % Plot the left- and right-running characteristics without data
     %   patch('Faces',faces_tri,'Vertices',nodes_tri/geom.yt,'FaceVertexCData',facesC_tri,...
     %         'EdgeColor','flat','FaceColor','none','LineWidth',1,'LineStyle','-');
     %   patch('Faces',faces_quad,'Vertices',nodes_quad/geom.yt,'FaceVertexCData',facesC_quad,...
     %         'EdgeColor','flat','FaceColor','none','LineWidth',1,'LineStyle','-');
     % end

%% NEW COLOR PLOTTING BY CLAUDE to get characteristic lines to appear. Lines commented out are from the original
     if (plot_colours)
         patch('Faces',faces_tri ,'Vertices',nodes_tri /geom.yt,'FaceVertexCData',nodesC_tri ,'FaceColor','interp','EdgeColor','k','LineWidth',0.25);
         patch('Faces',faces_quad,'Vertices',nodes_quad/geom.yt,'FaceVertexCData',nodesC_quad,'FaceColor','interp','EdgeColor','k','LineWidth',0.25);
         colormap(jet(256)); colorbar;
     else
         patch('Faces',faces_tri ,'Vertices',nodes_tri /geom.yt,'FaceVertexCData',nodesRGB_tri ,'EdgeColor','flat','FaceColor','none','LineWidth',1);
         patch('Faces',faces_quad,'Vertices',nodes_quad/geom.yt,'FaceVertexCData',nodesRGB_quad,'EdgeColor','flat','FaceColor','none','LineWidth',1);
     end

     % colormap(jet(256)) ;
     % colorbar; colorlim = clim;
     % 
     % % Do not plot the left- and right-running characteristics
     %   patch('Faces',faces_tri,'Vertices',nodes_tri/geom.yt,'FaceVertexCData',nodesC_tri,...
     %       'FaceColor','interp','LineStyle','none');
     % 
     %   patch('Faces',faces_quad,'Vertices',nodes_quad/geom.yt,'FaceVertexCData',nodesC_quad,...
     %         'FaceColor','interp','LineStyle','none');

  end

  axis equal;
  %axis([0 geom.xe 0 max(max(Y(1,:)))])
  %% Plot some arrows representing the flow direction at each point of intersection
  %quiver(X,Y,U,V)
  %saveas(fig1,'Characteristics.pdf')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot the static pressure on the wall + on the axis + comparison with 
  % 1D theoretical curve when the nozzle is chocked
  fig2=figure(2);
  condplot=zeros(size(X,1),size(X,2)); condplot = logical(condplot);
  wallnodesJ = geom.NI+1 : size(X,2) ;
  wallnodesJ = wallnodesJ ( X(1,wallnodesJ) < geom.xe ) ;
  for I=1:size(X,2) % Find the last index for each column -> this is the axis point
    axisnodesJ = find(U(:,I),1,'last');
    condplot( axisnodesJ , I ) = true;
  end
  
  [chocked] = get_Laval_theory(geom,params,X(1,wallnodesJ),Y(1,wallnodesJ)) ;
  semilogy(X(1,wallnodesJ)/geom.yt,pressure(1,wallnodesJ)/params.P,'r','linewidth',2); hold on; % Wall nodes
  semilogy(X(condplot)/geom.yt,pressure(condplot)/params.P,'b','linewidth',2); % Axis nodes
  semilogy(X(1,wallnodesJ)/geom.yt,chocked.pressure,'k','linewidth',2); % 1D theory on axis
  grid on; ylabel('Static pressure / Stagnation pressure [-]'); xlabel('Axial location x / y_t');
  legend('Wall','Axis','1D');
  axis([0 geom.xe 0.01 1]);
  set(gca,'xtick',0:1:10);
  fig = gcf;
  set(fig,'PaperUnits','normalized');
  set(fig,'PaperPosition',[0 0 1 0.4]);
  %saveas(fig2,'Pressure.pdf')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot the Mach number on the wall + on the axis + comparison with 
  % 1D theoretical curve when the nozzle is chocked
  fig3=figure(3);
  plot(X(1,wallnodesJ)/geom.yt,Mach(1,wallnodesJ),'r','linewidth',2); hold on; % Wall nodes
  plot(X(condplot)/geom.yt,Mach(condplot),'b','linewidth',2); % Axis nodes
  plot(X(1,wallnodesJ)/geom.yt,chocked.mach,'k','linewidth',2); % 1D theory on axis
  grid on; ylabel('Mach number [-]'); xlabel('Axial location x / y_t');
  legend('Wall','Axis','1D');
  axis([0 geom.xe 1 4]);
  set(gca,'xtick',0:1:10);
  fig = gcf;
  set(fig,'PaperUnits','normalized');
  set(fig,'PaperPosition',[0 0 1 0.4]);
  lgd=legend; set(legend,'Location','southeast');
  %saveas(fig3,'Mach_number.pdf')

end