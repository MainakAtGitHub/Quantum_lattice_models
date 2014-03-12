function [figure1,datarealmax]=plot_ldos_map_ft_path(ldosfile,scale,datarealmax,cut)
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
%figure1= figure('Position',[200, 50, 400, 300],'PaperUnits','centimeter','PaperPosition',[4 1 12 9]);
% Create axes
%axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1]);
% range=3
%axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'CameraViewAngle',3.88994451795861);
%range=5,'CameraViewAngle',1.6,
%axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',fntsz);
%grid(axes1,'on');
%hold(axes1,'all');
% do the fourier transform

kdividex=(max(xGridRange)-min(xGridRange))/RDiscrete(1);
kdividey=(max(yGridRange)-min(yGridRange))/RDiscrete(2);
x=X/RDiscrete(1);
y=Y/RDiscrete(2);
%xlim(axes1,[-2 2]);
%ylim(axes1,[-2 2]);
num=numel(localLdos);
% set up k-path
%[kx,ky]=meshgrid(-ceil(sqrt(RDiscrete(1)))*pi:pi/(kdividex+1):ceil(sqrt(RDiscrete(1)))*pi,-ceil(sqrt(RDiscrete(2)))*pi:pi/(kdividey+1):ceil(sqrt(RDiscrete(2)))*pi);
kx=-ceil(sqrt(RDiscrete(1)))*pi:pi/(kdividex+1):ceil(sqrt(RDiscrete(1)))*pi;
ky=kx;
szk=size(kx);
localLdosk=kx*0;
% fourier transform
for n=1:szk(1)
for m=1:szk(2)
localLdosk(n,m)=sum(sum(localLdos'.*exp(1i*(x*kx(n,m)+y*ky(n,m)))))/num;
end
end


plot(kx/pi,real(localLdosk));

if nargin < 3
    datarealmax=max(abs(localLdosk(:)));
else
    if isnan(datarealmax)
            datarealmax=max(abs(localLdosk(:)));
    end
end;
ylim([0, datarealmax]);
xlim([-2 2]);
%titleName = ['E = ', num2str(E), ', z = ', num2str(z), ' Bohr'];
%title(titleName);
xlabel('k_x/\pi=ky/\pi');
ylabel('FTdos' );
if isunix
    % create pdf of figure
    [~,filename,extension]=fileparts(ldosfile);
    print_pdf(['/tmp/',filename,extension,'ft_path.pdf']);
end;
