function p=plot_lattice_ldos_ft(ldosfile,plotN,datarealmax)
if nargin < 1
    ldosfile='lattice_greens_supercell_FeSe_N_15_M_9_U_0955_Vimp_5_E_minPt0084.mat';
end;
% input
scale='s';
fntsz=16;
lattice=false;
bluecolor=true;
fsz=20;
set(0,'DefaultAxesFontSize',fsz)
load(ldosfile,'-mat')
%lattice_greens_supercell_FeSe_N_15_M_9_U_0955_Vimp_5_E_minPt0084.mat
if ~(exist('N','var'))
    N = 11;
end;
if nargin <2
    plotN=N;
end;
xylabels=true;
if plotN<0
    plotN=-plotN;
    xylabels=false;
end;
if plotN>N
    plotN=N
end;
if ~(exist('E','var'))
    E = .0084;
end;
if ~(exist('nOrbitals','var'))
    nOrbitals = 10;
end;
if ~(exist('sublattice','var'))
    sublattice = 0;
end;


% total LDOS at Fe sites
latticeGreensDiag = diag(latticeGreens);
ldos = (-1/pi)*imag(reshape(latticeGreensDiag,nOrbitals*N,N));
Fe1LDOS = zeros(N,N);
Fe2LDOS = zeros(N,N);
if abs(sublattice)>0
for i = 1:N
    Fe1LDOS(i,:) = sum(ldos((i-1)*nOrbitals + (1:nOrbitals/2),:),1);
    Fe2LDOS(i,:) = sum(ldos((i-1)*nOrbitals + ((nOrbitals/2 + 1):nOrbitals),:),1);
end
Fe1LDOS = Fe1LDOS';
Fe2LDOS = Fe2LDOS';


% rotating coordinates by pi/4
n = ceil(N/2);
ldos2plot = diag(Fe1LDOS(:,n),0);
for i=1:(n-1)
    ldos2plot = ldos2plot + diag(Fe1LDOS((i+1):(N-i),n+i), 2*i) + diag(Fe1LDOS((i+1):(N-i),n-i), -2*i) + ....
    diag(Fe2LDOS(i:(N-i),n+i), 2*i-1) + diag(Fe2LDOS(i:(N-i),n-i+1), -(2*i-1));
end
else
    ldos2plot=ldos;
end;


num=numel(ldos2plot);
halfN=floor(N/2);
fine=2;
[kx,ky]=meshgrid(-pi:pi/(fine*halfN+1):pi,-pi:pi/(fine*halfN+1):pi);
szk=size(kx);
ldosk=ldos2plot*0;
range=-halfN:halfN;
[x,y]=meshgrid(range,range);
for n=1:szk(1)
    for m=1:szk(2)
        ldosk(n,m)=sum(sum(ldos2plot.*exp(1i*(x*kx(n,m)+y*ky(n,m)))))/num;
    end
end
figure1=figure;

if nargin < 3
    datarealmax=max(abs(ldosk(:)));
else
    if isnan(datarealmax)
            datarealmax=max(abs(ldosk(:)));
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
    %contourf(kx/pi,ky/pi,sqrt(real(ldosk)),'LineStyle','none');
    % Create surface
    % Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',16);
    box(axes1,'on');
hold(axes1,'all');
surface('Parent',axes1,'ZData',0*kx/pi,'YData',ky/pi,'XData',kx/pi,...
    'LineStyle','none',...
    'CData',sqrt(abs(real(ldosk))));
       % pcolor(kx/pi,ky/pi,sqrt(abs(real(ldosk))));%,'LineStyle','none');

    ticks=sign(tx).*sqrt(mtix*abs(tx));
    datarealmax=sqrt(datarealmax);
else
    %contourf(kx/pi,ky/pi,real(ldosk),'LineStyle','none');
        pcolor(kx/pi,ky/pi,real(ldosk));%,'LineStyle','none');

    datarealmax=datarealmax;
end;

%surfc(kx/pi,ky/pi,real(ldosk))
axis square
xlabel('k_x/\pi')
ylabel('k_y/\pi')
 bluemap(figure1)
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
%view(axes1,[0 90]);
caxis([-eps,datarealmax])
    set(h, 'YTick', ticksres1*datarealmax/256);
set(h, 'YTickLabel', labels1); 
if isunix
    % create pdf of figure
    [~,filename,extension]=fileparts(ldosfile);
   % if ~showcolorbar
    %                set(cb,'visible','off');
   % end;
    print_pdf(['/tmp/',filename,extension,'ft.pdf'])
    %print_eps(['/tmp/',filename,extension,'cut',num2str(N),'.eps'])
    
end;
