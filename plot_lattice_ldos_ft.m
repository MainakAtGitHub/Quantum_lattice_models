function figure1=plot_lattice_ldos_ft(ldosfile,scale,fine,fast,datarealmax)
if nargin < 1
    ldosfile='lattice_greens_supercell_FeSe_N_15_M_9_U_0955_Vimp_5_E_minPt0084.mat';
end;
% input
if nargin <2
    scale='';
end;
if nargin <3
    fine=1;
end;
if nargin <4
    fast=true;
end;
remove_bragg=true;
fntsz=16;
lattice=false;
bluecolor=true;
fsz=20;
set(0,'DefaultAxesFontSize',fsz)
load(ldosfile,'-mat')
%lattice_greens_supercell_FeSe_N_15_M_9_U_0955_Vimp_5_E_minPt0084.mat
if ~(exist('E','var'))
    E = .0084;
end;
[~, E,~]=getnumber(ldosfile,'_e_');
%E=-E
if ~(exist('nOrbitals','var'))
    nOrbitals = 10;
end;
if ~(exist('sublattice','var'))
    sublattice = 0
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

% remove the bragg peak!
if remove_bragg
    ldos2plot=ldos2plot-mean(ldos2plot(:));
end;
num=numel(ldos2plot);
halfN=floor(N/2);
kgrid=2*pi*(0:1/(fine*N):1-1/(fine*N));
[kx,ky]=meshgrid(kgrid,kgrid);
szk=size(kx);
ldosk=ldos2plot*0;
range=(1:N)-round(N/2);%-halfN:halfN;
[x,y]=meshgrid(range,range);
% Use fast fourier transform ?
% doesn't work 
%localLdos=wextend('2D','zpd',localLdos,(multiply-1)*szlLdos(1));
%ldos2plot=repmat(ldos2plot,fine,fine);
if ~fast
    % shift the position (doesn't matter if one takes the absolute value
    % later)
    % ldos2plot=fftshift(ldos2plot);
 for m=1:szk(1)
     for n=1:szk(2)
         ldosk(n,m)=sum(sum(ldos2plot.*exp(1i*(x*kx(n,m)+y*ky(n,m)))))/num/fine^2;
     end
 end
else
    %ldos2plot=repmat(ldos2plot,fine,fine);
    ldosk = fft2(ldos2plot,round(sqrt(num)*fine),round(sqrt(num)*fine))/num/fine^2;
   % ldosk = fftshift(ldosk);
end;
% do some cutoff of the k=0 component (not needed any more, see above)
if ~remove_bragg
ldosk(1,1)=0;
ldosk(1,1)=max(ldosk(:));
end;
ldosk = fftshift(ldosk);

figure1=figure('Position',[200, 50, 400, 300],'PaperUnits','centimeter','PaperPosition',[4 1 12 9]);
if nargin < 5
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
% extend the result to plot the (pi,pi) points etc. on both boundaries.
kgrid=2*pi*(0:1/(fine*N):1);
[kx1,ky1]=meshgrid(kgrid,kgrid);
% do also a shift for the k-matrices to get the labels right
kx=kx1-pi;
ky=ky1-pi;
% extension
ldosk=wextend('2d','ppd',ldosk,1);
% remove again the additional extension on the lower boundary
ldosk=ldosk(2:end,2:end);
% Create axes
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',16);
    box(axes1,'on');
hold(axes1,'all');
if scale=='s'
ldosk=sqrt(abs(ldosk));
    ticks=sign(tx).*sqrt(mtix*abs(tx));
    datarealmax=sqrt(datarealmax);
else
    %contourf(kx/pi,ky/pi,real(ldosk),'LineStyle','none');
      %  pcolor(kx/pi,ky/pi,abs(ldosk));%,'LineStyle','none');
ldosk=abs(ldosk);
    %datarealmax=datarealmax;
end;

surface('Parent',axes1,'ZData',0*kx/pi,'YData',ky/pi,'XData',kx/pi,...
    'LineStyle','none',...
    'CData',ldosk);

%surfc(kx/pi,ky/pi,real(ldosk))
axis square
xlabel('q_x/\pi')
ylabel('q_y/\pi')
%bma_map(figure1);
%blackmap(figure1);
 %bluemap(figure1)
 %hanaguri_map(figure1);
 hoffman_map(figure1);
 %fujita_map(figure1);
 axisshow=false
 if ~axisshow
     set(axes1, 'Visible','off')
 end;
 cbar=false;
 if cbar
     h = colorbar;
 end
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
if cbar
    set(h, 'YTick', ticksres1*datarealmax/256);
set(h, 'YTickLabel', labels1);
end;
if fine >1
    fn=['fine_',num2str(fine)];
else
    fn='';
end;
% % put some labels
%  annotation(figure1,'textbox',...
%      [0.2 0.8 0.4 0.1],...
%      'String',{['E=',num2str(E*1000),' meV']},...
%      'FontSize',20,...
%      'FitBoxToText','off',...
%      'EdgeColor','none', 'Color',[1 0 0]);
plotoctett=true;
if plotoctett
    symm=true;
[kx,ky]=banana_1band('~/itp/docs/real_space/BdG/bscco/Z3/input_SC_U015_N35_Z3.txt',E);
[o,figure1]=plot_octett_1band(kx,ky,E,figure1,symm);
else
%   % put some labels
%  annotation(figure1,'textbox',...
%      [0.2 0.8 0.4 0.1],...
%      'String',{['E=',num2str(E*1000),' meV']},...
%      'FontSize',20,...
%      'FitBoxToText','off',...
%      'EdgeColor','none', 'Color',[1 0 0]);  
end;
xlim([-2 2]);
ylim([-2 2]);
if isunix
    % create pdf of figure
    [~,filename,extension]=fileparts(ldosfile);
   % if ~showcolorbar
    %                set(cb,'visible','off');
   % end;
   if fine >2
       %print('-djpeg', ['/tmp/',filename,extension,fn,'_ft.jpg'],'-r200');
       print('-dpng', ['/tmp/',filename,extension,fn,'_ft.png'],'-r200');
   else
       print_pdf(['/tmp/',filename,extension,fn,'ft.pdf'])
   end;
    %print_eps(['/tmp/',filename,extension,'cut',num2str(N),'.eps'])
    
end;
plotcut=false;
if plotcut
% some plot of a cut for comparison
ldosk0 = fft2(ldos2plot);
if ~remove_bragg
ldosk0(1,1)=0;
ldosk0(1,1)=max(ldosk0(:));
end;
ldosk0 = fftshift(ldosk0)/num/fine^2;
if scale=='s'
ldosk0=sqrt(abs(ldosk0));
else
    %contourf(kx/pi,ky/pi,real(ldosk),'LineStyle','none');
      %  pcolor(kx/pi,ky/pi,abs(ldosk));%,'LineStyle','none');
ldosk0=abs(ldosk0);
    %datarealmax=datarealmax;
end;
sizek0=size(ldosk0);
sizek=size(kx);
kgrid0=2*pi*(1/(2*N):1/(N):1)-pi;


% Create figure
fig2 = figure;

% Create axes
axes1 = axes('Parent',fig2);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot((kx(fix(sizek(1)/2),1:end-1)+pi/(fine*N))/pi,ldosk(fix(sizek(1)/2),1:end-1),'Parent',axes1,'Marker','.','DisplayName','QPI (zero padding)');

% Create plot
plot(kgrid0/pi,ldosk0(ceil(sizek0(1)/2),:),'Parent',axes1,'MarkerFaceColor',[1 0 0],'Marker','square',...
    'DisplayName','QPI (bare)');

% Create legend
legend(axes1,'show');
end;


