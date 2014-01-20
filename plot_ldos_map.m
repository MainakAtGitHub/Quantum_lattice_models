function p=plot_ldos_map(ldosfile)
% to be modified for differen Wannier mesh
RDiscrete = [40 40 80];
load wannier_FeSe_4d_matrix_v2
load(ldosfile,'-mat')
zGridPoint = 0;
E = .0084;
localLdos = loacalLdos;
a = 7.23;
fntsz=16;
%z = zGrid(41 + zGridPoint);
k = strfind(ldosfile, '_z_');
ke = strfind(ldosfile, '_e_');
zstring1=ldosfile(k+3:length(ldosfile));
estring1=ldosfile(ke+3:k-1);
% works for both the normal calculation and the diag calculation
% should be programmed more sofisticated
try
    if zstring1(1)=='-'
        zpos=-str2num(num2str(zstring1(2:sum(isstrprop(zstring1, 'digit'))+1)));
    else
                zpos=str2num(num2str(zstring1(1:sum(isstrprop(zstring1, 'digit')))));
    end
catch
    zpos=input('could not find correct z-value, please enter: ','s');
end;
try
        E=str2num(estring1);
catch
    E=input('could not find correct energy value, please enter: ','s');
end;
[X, Y] = meshgrid(xGridRange,yGridRange);
% Create figure
figure1= figure('Position',[200, 50, 400, 300],'PaperUnits','centimeter','PaperPosition',[4 1 12 9]);
% Create axes
%axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1]);
% range=3
%axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'CameraViewAngle',3.88994451795861);
%range=5
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'CameraViewAngle',2.8,'FontSize',fntsz);
xlim(axes1,[min(X(:)) max(X(:))]);
ylim(axes1,[min(Y(:)) max(Y(:))]);
grid(axes1,'on');
hold(axes1,'all');
surf(X,Y,localLdos','LineStyle','none','FaceColor','interp');
% view from 1 Fe zone!
view(axes1,[45 90]);
 bluemap(figure1)
%view([0 90])
%pcolor(X,Y,localLdos');
%axis('square')
colorbar;
%titleName = ['E = ', num2str(E), ', z = ', num2str(z), ' Bohr'];
%title(titleName);
xlabel('x (Bohr)');
ylabel('y (Bohr)');
thickness=0.4;
ratio=0.8;
lgth=12;
offset=[-10,10];
x=[0,lgth]+offset(1);
y=[0,lgth]+offset(2);
z=[1 1]*max(localLdos(:));
h1=arrow3d(x,y,z,ratio,thickness);
text(1.1*x(2),1.1*y(2),z(2),'$x$','FontSize',1.5*fntsz,'Interpreter','latex')
set(h1,'facecolor',[1 0 0])
x=[0,-lgth]+offset(1);
y=[0,lgth]+offset(2);
z=[1 1]*max(localLdos(:));
h1=arrow3d(x,y,z,ratio,thickness);
text(1.1*x(2),1.1*y(2),z(2),'$y$','FontSize',1.5*fntsz,'Interpreter','latex')
% label the coordinate axes
set(h1,'facecolor',[1 0 0])
zposstring=['z=',sprintf('%1.3G',zpos/RDiscrete(3)),' c'];
% Create textbox
annotation(figure1,'textbox',...
        [0.0415264090747962 0.921052631578947 0.473511184910166 0.0690451293797412],...
    'String',{zposstring},...
    'FitBoxToText','off',...
    'LineStyle','none','FontSize',fntsz);

energystring=['\omega=',sprintf('%1.3G',E*1000),' meV'];
% Create textbox
annotation(figure1,'textbox',...
        [0.515068670435391 0.9375 0.353352382196187 0.0404822520324435],...
    'String',{energystring},...
    'FitBoxToText','off',...
    'LineStyle','none','FontSize',fntsz);
 set(gca,'FontSize', 16);
if isunix
    % create pdf of figure
    [~,filename,extension]=fileparts(ldosfile);
    print('-djpeg', ['/tmp/',filename,extension,'.jpg'],'-r200');
%print_pdf(['/tmp/',filename,extension,'cut',num2str(N),'.pdf'])
end;