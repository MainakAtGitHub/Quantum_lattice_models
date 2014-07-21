function [figure1,datarealmax]=plot_ldos_map(ldosfile,scale,datarealmax,cut,rotated)
if nargin < 2
    % default no sqrt scale!
    scale=''
end;
if nargin <4
    cut=5
end;
if nargin <5
    rotated='r';
end;
nolabel=false;
if cut<0
    nolabel=true;
    cut=-cut;
end;
% to be modified for different Wannier mesh
RDiscrete = [40 40 80]; % default value, correct value should be in "ldosfile"
%if ~(exist('wannier_filename','var'))
%    wannier_filename='wannier_FeSe_4d_matrix_v2.mat';
%        wannier_filename='./bscco/tom_input/WanF_Bi_2Sr_2CaCu_2O_8_vac_100_100_100_ReIm.out_conv.mat';
%end;
%load(wannier_filename);
load(ldosfile,'-mat')
RDiscrete
axistype='arrows';
axistype='lines';
zGridPoint = 0;
E = .0084;
if exist('loacalLdos','var')
    localLdos = loacalLdos;  
    clear loacalLdos;
end;
a = 7.23;
fntsz=16;
%z = zGrid(41 + zGridPoint);
% take out numbers from filename
[~, zpos, ~]=getnumber(ldosfile,'_z_');
[~, E,~]=getnumber(ldosfile,'_e_');

[X, Y] = meshgrid(xGridRange,yGridRange);
% Create figure
figure1= figure('Position',[200, 50, 400, 300],'PaperUnits','centimeter','PaperPosition',[4 1 12 9]);
% Create axes
%axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1]);
% range=3
%axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'CameraViewAngle',3.88994451795861);
%range=5
if ~isequal(rotated,'')
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'CameraViewAngle',1.6,'FontSize',fntsz);
else
    axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',fntsz);
end;
xlim(axes1,[min(X(:)) max(X(:))]);
ylim(axes1,[min(Y(:)) max(Y(:))]);
grid(axes1,'on');
hold(axes1,'all');
if nargin < 3 | isnan(datarealmax)
    datarealmax=max(abs(localLdos(:)));
else
    if isnan(datarealmax)
            datarealmax=max(abs(localLdos(:)));
    end
end;
lm=log(datarealmax)/log(10);
mtix=10^(ceil(lm));
% do some refinement to avoid only single labels
if (ceil(lm)-lm > 0.5)
    tx=[0:.025:.5]*2;
elseif    (ceil(lm)-lm > 0.35)
    tx=[0:0.05:1];
else
    tx=[0:0.1:1];
end;
ticks=mtix*tx; 
labels = num2str(repmat(sign(ticks).*(abs(ticks)), 1, 1)', 2);
if scale=='s'
    surf(X,Y,sqrt(abs(localLdos')),'LineStyle','none','FaceColor','flat');
    ticks=sign(tx).*sqrt(mtix*abs(tx));
    datarealmax=sqrt(datarealmax);
elseif scale =='l'
        surf(X,Y,log(abs(localLdos')),'LineStyle','none','FaceColor','flat');
            ticks=sign(tx).*log(mtix*abs(tx));
    datarealmax=log(datarealmax);
else
    surf(X,Y,abs(localLdos'),'LineStyle','none','FaceColor','flat');
    %datarealmax=datarealmax;
end;
global map
if isequal(map,'');
    bluemap(figure1);
elseif isequal(map,'neg');
    neg_bluemap(figure1);  
elseif isequal(map,'bma');
    bma_map(figure1);
elseif isequal(map,'song');
    song_map(figure1);
end;

%view([0 90])
%pcolor(X,Y,localLdos');
%axis('square')
    h = colorbar;
    ticks_res=round(ticks/datarealmax*256);
    % eliminate the same ticks_res
    ticksres1=ticks_res(1);
    labels1=labels(1,:);
    for n=2:length(ticks_res)
        if ticks_res(n)> ticksres1(length(ticksres1))
            ticksres1=[ticksres1,ticks_res(n)];
            labels1=[labels1;labels(n,:)];
        end
    end;
    if ~isequal(scale,'l')
        caxis([0,datarealmax])
        allAxesInFigure = findall(figure1,'type','axes');
        set(allAxesInFigure,'CLim',[0 datarealmax],'FontSize',fntsz); 
    end;
% view from 1 Fe zone!
if rotated=='r'
    view(axes1,[45 90]);
elseif rotated=='n'
    view([0 90])
elseif rotated==''
    axistype='';
end;
redlines=true;
    if ~isequal(scale,'l')

caxis([-eps,datarealmax])
    end
set(h, 'YTick', ticksres1*datarealmax/256);
set(h, 'YTickLabel', labels1);
%titleName = ['E = ', num2str(E), ', z = ', num2str(z), ' Bohr'];
%title(titleName);
%xlabel('x (Bohr)');
%ylabel('y (Bohr)');
xlabel('x');
ylabel('y');
thickness=0.4;
ratio=0.8;
lgth=0.3*RDiscrete(1);
if ~exist('sublattice','var')
    sublattice=1
end;
if sublattice==1
    offset(1)=-0.25*RDiscrete(1);
    offset(2)=0.25*RDiscrete(2);
else
    offset=[0 0];
end;
z=[1 1]*datarealmax;
switch axistype
    case 'arrows'
        % insert coordinate system manually (arrows)
        x=[0,lgth]+offset(1);
        y=[0,lgth]+offset(2);
        z=[1 1]*datarealmax;
        h1=arrow3d(x,y,z,ratio,thickness);
        text(1.1*x(2),1.1*y(2),z(2),'$x$','FontSize',1.5*fntsz,'Interpreter','latex')
        set(h1,'facecolor',[1 0 0])
        x=[0,-lgth]+offset(1);
        y=[0,lgth]+offset(2);
        h1=arrow3d(x,y,z,ratio,thickness);
        text(1.1*x(2),1.1*y(2),z(2),'$y$','FontSize',1.5*fntsz,'Interpreter','latex')
        % put in z-component and energy as text
        set(h1,'facecolor',[1 0 0])
    case 'lines'
        switch rotated
            case 'r'
                x=cut/4*[-RDiscrete(1),RDiscrete(2)]+offset(1)+.5;
                y=cut/4*[-RDiscrete(1),RDiscrete(2)]+offset(2)+.5;
                if redlines
                    plot3(x,y,z,'r');
                end;
                x=cut/4*[RDiscrete(1),-RDiscrete(2)]+offset(1)+.5;
                y=cut/4*[-RDiscrete(1),RDiscrete(2)]+offset(2)+.5;
                if redlines
                    plot3(x,y,z,'r');
                end;
                if exist('cut','var')
                    % plot black box with corresponding cut and cuttof the image with white
                    % boxes
                    boxx=[-1, 0, 1, 0,-1];
                    boxy=[ 0, 1, 0,-1, 0];
                    z=[1 1 1 1 1]*datarealmax*1.02;
                    boxx=cut/2*RDiscrete(1)*boxx+offset(1);
                    boxy=cut/2*RDiscrete(2)*boxy+offset(2);
                    plot3(boxx,boxy,z*1.02,'k');
                    boxscale=3;
                    boxcolor='white';
                    % plot white box
                    boxx1=[ 0 1 boxscale  0 0]*cut/2*RDiscrete(1)+offset(1);
                    boxy1=[-1 0 0 -boxscale -1]*cut/2*RDiscrete(2)+offset(2);
                    patch(boxx1,boxy1,z,boxcolor,'LineStyle','none');
                    boxx1=[-1  0  0 -boxscale -1]*cut/2*RDiscrete(1)+offset(1);
                    boxy1=[ 0 -1 -boxscale  0  0]*cut/2*RDiscrete(2)+offset(2);
                    patch(boxx1,boxy1,z,boxcolor,'LineStyle','none');
                    boxx1=[0 -1 -boxscale 0 0]*cut/2*RDiscrete(1)+offset(1);
                    boxy1=[1  0  0 boxscale 1]*cut/2*RDiscrete(2)+offset(2);
                    patch(boxx1,boxy1,z,boxcolor,'LineStyle','none');
                    boxx1=[1 0 0 boxscale 1]*cut/2*RDiscrete(1)+offset(1);
                    boxy1=[0 1 boxscale 0 0]*cut/2*RDiscrete(2)+offset(2);
                    patch(boxx1,boxy1,z,boxcolor,'LineStyle','none');
                    xlim(boxscale*RDiscrete(1)*cut/2*[ -1 1]);
                    ylim(boxscale*RDiscrete(2)*cut/2*[ -1 1]);
                end
            case 'n'
                x=cut/2*[-RDiscrete(1),RDiscrete(1)]+offset(1);
                y=cut/2*[0,0]+offset(2);
                %plot3(x,y,z,'r');
                x=cut/2*[0,0]+offset(1);
                y=cut/2*[-RDiscrete(2),RDiscrete(2)]+offset(2);
               % plot3(x,y,z,'r');
                if exist('cut','var')
                    % plot black box with corresponding cut and cuttof the image with white
                    % boxes
                    boxx=[-1, 1, 1, -1,-1];
                    boxy=[ -1, -1, 1,1, -1];
                    z=[1 1 1 1 1]*datarealmax*1.02;
                    boxx=cut/2*RDiscrete(1)*boxx+offset(1);
                    boxy=cut/2*RDiscrete(2)*boxy+offset(2);
                    plot3(boxx,boxy,z*1.02,'k');
                    boxscale=3;
                    boxcolor='white';
                    % plot white box
                    boxx1=[ -1 1 boxscale  -boxscale -1]*cut/2*RDiscrete(1)+offset(1);
                    boxy1=[-1 -1 -boxscale -boxscale -1]*cut/2*RDiscrete(2)+offset(2);
                    patch(boxx1,boxy1,z,boxcolor,'LineStyle','none');
                    boxx1=[-1  -1  -boxscale -boxscale -1]*cut/2*RDiscrete(1)+offset(1);
                    boxy1=[ -1 1 boxscale  -boxscale  -1]*cut/2*RDiscrete(2)+offset(2);
                    patch(boxx1,boxy1,z,boxcolor,'LineStyle','none');
                    boxx1=[-1 1 boxscale -boxscale -1]*cut/2*RDiscrete(1)+offset(1);
                    boxy1=[1  1  boxscale boxscale 1]*cut/2*RDiscrete(2)+offset(2);
                    patch(boxx1,boxy1,z,boxcolor,'LineStyle','none');
                    boxx1=[1 1 boxscale boxscale 1]*cut/2*RDiscrete(1)+offset(1);
                    boxy1=[1 -1 -boxscale boxscale 1]*cut/2*RDiscrete(2)+offset(2);
                    patch(boxx1,boxy1,z,boxcolor,'LineStyle','none');
                    xlim(boxscale*RDiscrete(1)*cut/sqrt(2)*[ -1 1]);
                    ylim(boxscale*RDiscrete(2)*cut/sqrt(2)*[ -1 1]);
                end;
                
        end;
end
if ~nolabel
zposstring=['z=',sprintf('%1.3G',zpos/RDiscrete(3)),' c'];
zposstring='';
% Create textbox
annotation(figure1,'textbox',...
        [0.08 0.921052631578947 0.473511184910166 0.0690451293797412],...
    'String',{zposstring},...
    'FitBoxToText','off',...
    'LineStyle','none','FontSize',fntsz);

energystring=['\omega=',sprintf('%1.3G',E*1000),' meV'];
energystring='';
% Create textbox
annotation(figure1,'textbox',...
        [0.515068670435391 0.9375 0.353352382196187 0.0404822520324435],...
    'String',{energystring},...
    'FitBoxToText','off',...
    'LineStyle','none','FontSize',fntsz);
 set(gca,'FontSize', 16);
end;
if isunix
    % create pdf of figure
    [~,filename,extension]=fileparts(ldosfile);
    if nolabel
            set(h,'visible','off');
            print('-dpng', ['/tmp/',filename,extension,'.png'],'-r150');
    else
set(gcf, 'Renderer', 'OpenGL');
%set(gcf, 'Renderer', 'zbuffer');
print('-djpeg', ['/tmp/',filename,extension,'.jpg'],'-r150');
print('-dpng', ['/tmp/',filename,extension,'.png'],'-r150');

    end;
%print_pdf(['/tmp/',filename,extension,'cut',num2str(N),'.pdf'])
if nolabel
    % print also just the colorbar
    childr = get(axes1,'children');
    set(childr,'visible','off');
    set(axes1,'visible','off');
    set(h,'visible','on');
    %print('-djpeg', ['/tmp/',filename,extension,'colorbar.jpg'],'-r200');
    print_pdf(['/tmp/',filename,extension,'colorbar.pdf']);
    set(childr,'visible','on');
    set(axes1,'visible','on');
    set(h,'visible','on');
end
end;