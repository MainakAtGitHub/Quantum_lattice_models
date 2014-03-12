function [figure1,datarealmax]=plot_ldos_map_ft(ldosfile,scale,datarealmax,cut)
if nargin < 2
    % default no sqrt scale!
    scale=''
end;
if nargin <4
    cut=5
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
%range=5,'CameraViewAngle',1.6,
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',fntsz);
grid(axes1,'on');
hold(axes1,'all');
% do the fourier transform

kdividex=(max(xGridRange)-min(xGridRange))/RDiscrete(1);
kdividey=(max(yGridRange)-min(yGridRange))/RDiscrete(2);
x=X/RDiscrete(1);
y=Y/RDiscrete(2);
xlim(axes1,[-2 2]);
ylim(axes1,[-2 2]);
num=numel(localLdos);
[kx,ky]=meshgrid(-ceil(sqrt(RDiscrete(1)))*pi:pi/(kdividex+1):ceil(sqrt(RDiscrete(1)))*pi,-ceil(sqrt(RDiscrete(2)))*pi:pi/(kdividey+1):ceil(sqrt(RDiscrete(2)))*pi);
szk=size(kx);
localLdosk=kx*0;
for n=1:szk(1)
for m=1:szk(2)
localLdosk(n,m)=sum(sum(localLdos'.*exp(1i*(x*kx(n,m)+y*ky(n,m)))))/num;
end
end




if nargin < 3
    datarealmax=max(abs(localLdosk(:)));
else
    if isnan(datarealmax)
            datarealmax=max(abs(localLdosk(:)));
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
    surf(kx/pi,ky/pi,sqrt(abs(localLdosk)),'LineStyle','none','FaceColor','interp');
    ticks=sign(tx).*sqrt(mtix*abs(tx));
    datarealmax=sqrt(datarealmax);
else
    surf(kx/pi,ky/pi,abs(localLdosk),'LineStyle','none','FaceColor','interp');
    datarealmax=datarealmax;
end;
 bluemap(figure1)


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
 caxis([0,datarealmax])
        allAxesInFigure = findall(figure1,'type','axes');
        set(allAxesInFigure,'CLim',[0 datarealmax],'FontSize',fntsz); 
view(axes1,[0 90]);
caxis([-eps,datarealmax])
    set(h, 'YTick', ticksres1*datarealmax/256);
set(h, 'YTickLabel', labels1);
%titleName = ['E = ', num2str(E), ', z = ', num2str(z), ' Bohr'];
%title(titleName);
xlabel('k_x/\pi');
ylabel('k_y/\pi' );
if isunix
    % create pdf of figure
    [~,filename,extension]=fileparts(ldosfile);
    if nolabel
            set(h,'visible','off');
            print('-dpng', ['/tmp/',filename,extension,'_ft.png'],'-r200');
    else
        print('-djpeg', ['/tmp/',filename,extension,'_ft.jpg'],'-r200');
    end;
%print_pdf(['/tmp/',filename,extension,'cut',num2str(N),'.pdf'])
if nolabel
    % print also just the colorbar
    childr = get(axes1,'children');
    set(childr,'visible','off');
    set(axes1,'visible','off');
    set(h,'visible','on');
    %print('-djpeg', ['/tmp/',filename,extension,'colorbar.jpg'],'-r200');
    print_pdf(['/tmp/',filename,extension,'colorbar_ft.pdf']);
    set(childr,'visible','on');
    set(axes1,'visible','on');
    set(h,'visible','on');
end
end;
