function [figure1,datarealmax]=plot_ldos_map_ft(ldosfile,scale,fine,cut,datarealmax)
if nargin < 2
    % default no sqrt scale!
    scale=''
end;
if nargin <4
    cut=50
end;
if nargin < 3
    fine=1;
end;
nolabel=false;
if cut<0
    nolabel=true;
    cut=-cut;
end;
remove_bragg=true;

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
%E=-E
% cut the ldosmap to get rid of the bug!
cutoff=50
range1=(cutoff+1):(length(xGridRange)-cutoff-1);
localLdos=localLdos(range1,range1);
xGridRange=xGridRange(range1);
yGridRange=yGridRange(range1);
szlLdos=size(localLdos);
% doesn't work 
%fine=10;
%[X, Y] = meshgrid(xGridRange,yGridRange);
% Create figure
figure1= figure('Position',[200, 50, 400, 300],'PaperUnits','centimeter','PaperPosition',[4 1 12 9]);
% Create axes
%axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1]);
% range=3
%axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'CameraViewAngle',3.88994451795861);
%range=5,'CameraViewAngle',1.6,
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',fntsz);
%grid(axes1,'on');
box(axes1,'on');
hold(axes1,'all');
% do the fourier transform

kdividex=fine*(max(xGridRange)-min(xGridRange)+1)/RDiscrete(1)-1;
kdividey=fine*(max(yGridRange)-min(yGridRange)+1)/RDiscrete(2)-1;
%x=X/RDiscrete(1);
%y=Y/RDiscrete(2);
num=numel(localLdos);
%kgrid=2*pi*(0:1/(fine*N):1-1/(fine*N));
%[kx,ky]=meshgrid(kgrid,kgrid);
[kx,ky]=meshgrid(-ceil(sqrt(RDiscrete(1)))*pi:2*pi/(kdividex+1):ceil(sqrt(RDiscrete(1)))*pi,-ceil(sqrt(RDiscrete(2)))*pi:2*pi/(kdividey+1):ceil(sqrt(RDiscrete(2)))*pi);
szk=size(kx);
localLdosk=kx*0;
% for n=1:szk(1)
% for m=1:szk(2)
% localLdosk(n,m)=sum(sum(localLdos'.*exp(1i*(x*kx(n,m)+y*ky(n,m)))))/num;
% end
% end
% remove outer edges from input (already done above)
%localLdos=localLdos(1:end-1,1:end-1);
% to do: use the fast fourier transform from Matlab for speedup
if remove_bragg
    localLdos=localLdos-mean(localLdos(:));
    %remove also secondary peaks
    peakorder=2;
    while peakorder > 1
        peakorder=peakorder-1
        % do the FFT with the small lattice
        localLdosk1 = fft2(localLdos);
        % remove peaks at +/-2pi (4 largest values)
        [sortedValues,~] = sort(localLdosk1(:),'descend');  %# Sort the values in
        maxValues = sortedValues(1:4);  %# Get the 5 largest values
        maxIndex = ismember(localLdosk1,maxValues);     %# Get a logical index of all values
        %#   equal to the 5 largest values                                   %#   descending order
        % manual removal for special map, to
        % be adapted!
        if false
            localLdosk2=localLdosk1*0;
            localLdosk2(31)=localLdosk1(31);
            localLdosk2(571)=localLdosk1(571);
            localLdosk2(18001)=localLdosk1(18001);
            localLdosk2(342001)=localLdosk1(342001);
            localLdos=localLdos-ifft2(localLdosk2);
        else
            localLdosk1=localLdosk1.*maxIndex;
            localLdos=localLdos-ifft2(localLdosk1);
        end;
        %        localLdosk1=(max(abs(localLdosk1(:))*0.95)<abs(localLdosk1)).*localLdosk1;
    end,
end;
localLdosk = fft2(localLdos,fine*sqrt(num),fine*sqrt(num));
localLdosk = fftshift(localLdosk);
szldos=size(localLdosk);


localLdosk=wextend('2d','ppd',localLdosk,1);
% remove again the additional extension on the lower boundary
localLdosk=localLdosk(2:end,2:end);
if ~remove_bragg
% do some cutoff of the k=0 component
localLdosk(size(localLdosk,1)/2+0.5,size(localLdosk,2)/2+0.5)=0;
localLdosk(size(localLdosk,1)/2+0.5,size(localLdosk,2)/2+0.5)=max(localLdosk(:));
end;

xrang=(1:szk(1))+ (szldos(1)-szk(1)-1)/2+1;
yrang=(1:szk(2))+ (szldos(2)-szk(2)-1)/2+1;
localLdosk=localLdosk(xrang,yrang);

if nargin < 5
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
 %bluemap(figure1)
%bma_map(figure1);
%blackmap(figure1);
%hanaguri_map(figure1);
hoffman_map(figure1);
%fujita_map(figure1);

%view([0 90])
%pcolor(X,Y,localLdos');
%axis('square')
 axisshow=false
 if ~axisshow
     set(axes1, 'Visible','off')
 end;
cbar=false
if cbar
    h = colorbar;
end;
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
if cbar
    set(h, 'YTick', ticksres1*datarealmax/256);
set(h, 'YTickLabel', labels1);
end;
%titleName = ['E = ', num2str(E), ', z = ', num2str(z), ' Bohr'];
%title(titleName);
xlabel('q_x/\pi');
ylabel('q_y/\pi' );
xlim(axes1,[-2 2]);
ylim(axes1,[-2 2]);
if fine >1
    fn=['fine_',num2str(fine)];
else
    fn='';
end;
shading flat

plotoctett=false;
if plotoctett
    symm=true;
[kx,ky]=banana_1band('~/itp/docs/real_space/BdG/bscco/Z3/input_SC_U015_N35_Z3.txt',E);
h=datarealmax;
[o,figure1]=plot_octett_1band(kx,ky,E,figure1,symm,h);
else
  % put some labels
%  annotation(figure1,'textbox',...
%      [0.2 0.8 0.4 0.1],...
%      'String',{['E=',num2str(E*1000),' meV']},...
%      'FontSize',20,...
%      'FitBoxToText','off',...
%      'EdgeColor','none', 'Color',[1 0 0]);  
end;
if isunix
    % create pdf of figure
    [~,filename,extension]=fileparts(ldosfile);
    if nolabel
        set(h,'visible','off');
        print('-dpng', ['/tmp/',filename,extension,fn,'_ft.png'],'-r200');
        xlim(axes1,[-1 1]);
        ylim(axes1,[-1 1]);
        print('-dpng', ['/tmp/',filename,extension,fn,'_ft_zoom.png'],'-r200');
        xlim(axes1,[-2 2]);
        ylim(axes1,[-2 2])
    else
        %print('-djpeg', ['/tmp/',filename,extension,fn,'_ft.jpg'],'-r200');
        print('-dpng', ['/tmp/',filename,extension,fn,'_ft.png'],'-r200');
        %xlim(axes1,[-1 1]);
        %ylim(axes1,[-1 1]);
        %print('-djpeg', ['/tmp/',filename,extension,fn,'_ft_zoom.jpg'],'-r200');
        xlim(axes1,[-2 2]);
        ylim(axes1,[-2 2])        
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
