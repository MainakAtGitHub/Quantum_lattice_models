function p=plot_lattice_ldos(ldosfile,plotN,datarealmax)
if nargin < 1
    ldosfile='lattice_greens_supercell_FeSe_N_15_M_9_U_0955_Vimp_5_E_minPt0084.mat';
end;
% input
lattice=false;
bluecolor=true;
fsz=20;
set(0,'DefaultAxesFontSize',fsz)
load(ldosfile,'-mat')
%lattice_greens_supercell_FeSe_N_15_M_9_U_0955_Vimp_5_E_minPt0084.mat
if ~(exist('N','var'))
    N = 11;
end;
xylabels=true;
if nargin <2
    plotN=N;
    xylabels=false;
end;
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

figure1=figure;
axes1 = axes('Parent',figure1,'YDir','reverse',...
    'PlotBoxAspectRatio',[1 1 1],...
    'Layer','top');
 %   'XTickLabel',{'-5','-4','-3','-2','-1','0','1','2','3','4','5'},...
  %  'YTickLabel',{'-5','-4','-3','-2','-1','0','1','2','3','4','5'},...
%hold(axes1,'all');
%box(axes1,'on');

pcolorplot=false;
cptn='LDOS [1/eV]';
% zoom the figure to show only plotN points
diffN=N-plotN;
if sublattice==1
    ldos2plot=ldos2plot(1+diffN/2:N-diffN/2,1+diffN/2:N-diffN/2);
elseif sublattice==-1
    % not correct yet (to be done)
    ldos2plot=ldos2plot(-1+diffN/2:N-diffN/2,-1+diffN/2:N-diffN/2);
else
    ldos2plot=ldos2plot(diffN/2+1:N-diffN/2,diffN/2+1:N-diffN/2);
end
N=plotN;
n = ceil(N/2);
showcolorbar=true;
if nargin < 3
    datarealmax=max(ldos2plot(:));
end;
if datarealmax<0
    datarealmax=-datarealmax;
    showcolorbar=false;
end;
% plotting (using pcolor)
if pcolorplot
    ldos2plot = flipud(ldos2plot);
    ldos2plot(:,N+1) = 0;
    ldos2plot(N+1,:) = 0;
    pcolor(ldos2plot,'Parent',axes1);
else
    image(ldos2plot,'Parent',axes1,'CDataMapping','scaled');
    tickx_num=(-n+1:1:n-1);
    tickx=cellstr(num2str(tickx_num(:)));
    %tickx={'-5','-4','-3','-2','-1','0','1','2','3','4','5'};
    if ~xylabels
        for i=1:numel(tickx)
            tickx{i}='';
        end;
    end;
        label_boxes_ldos(numel(tickx),tickx);
        % move the labels out of the ticks
        dp=0.05*N;
        yh=get(axes1,'ylabel');
        posy=get(yh,'position');
        set(yh,'position',[posy(1)-dp posy(2)])
        yh=get(axes1,'xlabel');
        posy=get(yh,'position');
        set(yh,'position',[posy(1) posy(2)+dp])
        if xylabels
            xlabel('\Delta x')
            ylabel('\Delta y')
        end;      
end;
axis('square');
title(['\omega = ',num2str(E*1000), ' meV']);
% colorbar schemes

cb=colorbar;
if bluecolor
 bluemap(figure1)
end;
caxis([0,datarealmax])
caxis
% set caption to colorbar
if ~strcmp(cptn,'')
    zlab = get(cb,'ylabel');
    set(zlab,'String',cptn,'FontSize',fsz);
end;
if lattice
    [x,y]=meshgrid(1:N);
    hold on;
    pointsize=160;
    lnwth=0.6;
    % cut of the central point
    floor(N^2/2)
    x1=[x(1:floor(N^2/2)),x(floor(N^2/2)+2:N^2)];
    y1=[y(1:floor(N^2/2)),y(floor(N^2/2)+2:N^2)]; 
    scatter(x1(:),y1(:),pointsize,'MarkerEdgeColor','k',...
              'MarkerFaceColor','r',...
              'LineWidth',lnwth);
              [x,y]=meshgrid(1:N-1);
          xse=x(mod(x(:)+y(:),2)==0)+0.5;
          yse=y(mod(x(:)+y(:),2)==0)+0.5;
              scatter(xse(:),yse(:),pointsize,'^','MarkerEdgeColor','k',...
              'MarkerFaceColor','y',...
              'LineWidth',lnwth);
                    xse=x(mod(x(:)+y(:),2)==1)+0.5;
          yse=y(mod(x(:)+y(:),2)==1)+0.5;
              scatter(xse(:),yse(:),pointsize,'v','MarkerEdgeColor','k',...
              'MarkerFaceColor','y',...
              'LineWidth',lnwth);
          scatter(n ,n,pointsize*3.2,'h','MarkerEdgeColor','k',...
              'MarkerFaceColor','g',...
              'LineWidth',lnwth);
end;

if isunix
    % create pdf of figure
    [~,filename,extension]=fileparts(ldosfile);
    if ~showcolorbar
                    set(cb,'visible','off');
    end;
    print_pdf(['/tmp/',filename,extension,'cut',num2str(N),'.pdf'])
    %print_eps(['/tmp/',filename,extension,'cut',num2str(N),'.eps'])
    
end;
