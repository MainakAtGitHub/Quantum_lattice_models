function figure1=plot_FeSe_elementary(z)



% to be modified for differen Wannier mesh
RDiscrete = [40 40 80];
figure1= figure('Position',[200, 50, 150, 150]);
M=diag([7.13,7.13,10.44]);
showaxes=false;
% Create axes
if ~showaxes
    axes1 = axes('PlotBoxAspectRatio',[1 1 1], 'Visible','off','Parent',figure1);
else
    axes1 = axes('PlotBoxAspectRatio',[1 1 1],'Parent',figure1);
    view(axes1,[-14 18]);
    % does only work with diagonal M!
    xlim([0, M(1,1)])
        ylim([0, M(2,2)])
   zlim(0.5*[-M(3,3), M(3,3)])
   xlabel('x bohr')
   ylabel('y bohr');
   zlabel('z bohr');
end;

view(axes1,[0 0]);
hold(axes1,'all');
    
% from Tom's file Fe and Se positions
position= [0.25, 0.75, 0;0.75, 0.25, 0;0.75, 0.75, 0.265;0.25, 0.25, 0.735-1];
% M is diagonal, thus matrix multiplication simple
M=diag([7.13,7.13,10.44]);
positionx=position(:,1)*M(1,1);
positiony=position(:,2)*M(2,2);
positionz=position(:,3)*M(3,3);
% plot elementary cell
a1=[0 0 0 0 0]-0.5;
a1a=[1 1 1 1 1]-0.5;
a2=[0 1 1 0 0];
a3=[0 0 1 1 0];
lwthcube=2;
plot3(a2*M(1,1),a3*M(2,2),a1*M(3,3),'-k','Linewidth',lwthcube,'Marker','o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
plot3(a2*M(1,1),a3*M(2,2),a1a*M(3,3),'-k','Linewidth',lwthcube,'Marker','o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
b1=[0 1]-0.5;
one=[ 1 1];
plot3(one*0+eps,one*0+eps,b1*M(3,3),'-k','Linewidth',lwthcube);
plot3(one*M(1,1),one*0+eps,b1*M(3,3),'-k','Linewidth',lwthcube);
plot3(one*M(1,1),one*M(2,2),b1*M(3,3),'-k','Linewidth',lwthcube);
plot3(one*0+eps,one*M(2,2),b1*M(3,3),'-k','Linewidth',lwthcube);

% plot atoms
% plot line indicating cut position
    hold on;
    pointsize=80;
    lnwth=0.6;
    scatter3(positionx(1:2),positiony(1:2),positionz(1:2),pointsize,'MarkerEdgeColor','k',...
              'MarkerFaceColor','r',...
              'LineWidth',lnwth);
              scatter3(positionx(3),positiony(3),positionz(3),pointsize,'^','MarkerEdgeColor','k',...
              'MarkerFaceColor','y',...
              'LineWidth',lnwth);
              scatter3(positionx(4),positiony(4),positionz(4),pointsize,'v','MarkerEdgeColor','k',...
              'MarkerFaceColor','y',...
              'LineWidth',lnwth);
             % cut of the central point
             [xmesh,ymesh]=meshgrid(0:1,0:1);
             zmesh=ones(2,2)*M(3,3)*z/RDiscrete(3);
             surfc(xmesh*M(1,1),ymesh*M(2,2),zmesh,'LineWidth',3,'EdgeColor',[0 1 0]);
             daspect([1 1 1])
   zstring=sprintf('%02d',z+RDiscrete(3)/2)        
 print('-djpeg', ['/tmp/FeSe_elementary',zstring,'.jpg']);
 print_pdf(['/tmp/FeSe_elementary',zstring,'.pdf']);
             
