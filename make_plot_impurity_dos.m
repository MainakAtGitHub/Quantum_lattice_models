function [figure5,energy,sumorbitalLDOS]=make_plot_impurity_dos(inputfile,smoothenergy,plrange,plotall,omega0,figure5)
global plotrange homogeneous
if isempty(homogeneous)
    homogeneous=false
end;
if nargin <1
    inputfile='LDOS_FeSe_Milan_Gamma_Vimp_4_N_9_M_40_ita_0.003'
end;
% number of ldos positions to be plotted
if nargin <3
    plrange=inf;
end;
if nargin <4
    plotall=true;
end
if plrange<0
    plotall=false
end;
if nargin < 5
    omega0=[];
end
%plotall=true;
coloruf1=[250 	70 	22 ]/255;
% pure red instead of orange
coloruf1=[250 	0 	0 ]/255;
% pure blue instead of UF_blue
coloruf2= [0 	48 	135]/255;
% support for wildcards
kstar = strfind(inputfile, '*');
kquestion=strfind(inputfile, '?');
if sum(kstar,kquestion)>0
   inputfilelist=dir(inputfile);
   [directory,~,~]=fileparts(inputfile);
   if ~isempty(directory)
       directory=[directory,'/'];
   end;
   sza=size(inputfilelist);
   for numfile=1:sza(1);
       disp(['Proccessing ',directory,inputfilelist(numfile).name]);
       if nargin < 2
           make_plot_impurity_dos([directory,inputfilelist(numfile).name]);
       else
           make_plot_impurity_dos([directory,inputfilelist(numfile).name],smoothenergy);
       end
  
       if sza(1)>5
           close all;
       end;
   end;
   return;
end;
if nargin < 2
    if isempty(strfind(inputfile, 'tetra'))
        smoothenergy=0;
    else
        %default: smoothing of 2 meV if "tetra" in filename
        smoothenergy=0.002;
    end
end

if (isempty(plotrange))
    plotrange=[-0.6 0.6];
        %plotrange=[-2 2];
end;


fsz=14;
% load input
try
    load(inputfile,'-mat')
catch
    load(inputfile)
end;
% convert input to new format if needed
if ~exist('orbitalLDOS','var')
    szdata=size(orbitalLDOSFarAway);
    efforb=szdata(1);
    orbitalLDOS((1:efforb),:)=orbitalLDOSFarAway;
    orbitalLDOS(efforb+(1:efforb),:)=orbitalLDOSImp;
    orbitalLDOS(2*efforb+(1:efforb),:)=orbitalLDOSImpNN;
    orbitalLDOS(3*efforb+(1:efforb),:)=orbitalLDOSImpNNN;
    LDOSsites=[-inf,-inf; 0,0; 0,1;1,1];
end;
nDosSites = size(LDOSsites,1);
if homogeneous
    nDosSites=1
end;
orb={'orbital1','orbital2','orbital3','orbital4','orbital5','total'};
%orb={'d_{z^2}','d_{x^2-y^2}','d_{yz}','d_{xz}','d_{xy}','total'}; % labels for Tom's FeSe model
% fix for non existing variable
if ~exist('efforb','var')
    efforb=1;
end;
orb={orb{1:efforb},orb{6}};
% some smoothing if necessary
de=energy(2)-energy(1);
smooth=floor(smoothenergy/de);
    egrid=-(400*de):de:(400*de);
    %global temperature
    temperature=smoothenergy;%2/11400;
    dF=-fermi_prime_func(egrid/temperature);
    dF=dF/sum(dF);
for n=1:nDosSites
    %tmp=sg_smooth(orbitalLDOS((n-1)*efforb+(1:efforb),:),smooth);
    tmp=orbitalLDOS((n-1)*efforb+(1:efforb),:);
    % do some smoothing with temperature
    if temperature >0
    for sz=1:size(tmp,1)
        % repeat points on the boundaries
        tmp1=[repmat(tmp(sz,1),1,floor(numel(dF)/2)),tmp(sz,:),repmat(tmp(sz,end),1,floor(numel(dF)/2))];
    tmp(sz,:)=conv(tmp1',dF,'valid')';
    end
    tmp(tmp<0)=0;
    end
    orbitalLDOS((n-1)*efforb+(1:efforb),:)=tmp;
end;
clear tmp;
% orbitalLDOSFarAway=sg_smooth(orbitalLDOSFarAway,smooth);
% orbitalLDOSImp=sg_smooth(orbitalLDOSImp,smooth);
% orbitalLDOSImpNN=sg_smooth(orbitalLDOSImpNN,smooth);
% orbitalLDOSImpNNN=sg_smooth(orbitalLDOSImpNNN,smooth);
% % fix negative values (numerical error, smoothing artefacts)
% orbitalLDOSFarAway(orbitalLDOSFarAway<0)=0;
% orbitalLDOSImp(orbitalLDOSImp<0)=0;
% orbitalLDOSImpNN(orbitalLDOSImpNN<0)=0;
% orbitalLDOSImpNNN(orbitalLDOSImpNNN<0)=0;
% DOS far away (without impurity)

% get filename of inputfile without path (doesn't work in windows yet)
k = findstr(inputfile, '/');
inputfile_orig=inputfile;
if ~isempty(k)
    	inputfile=inputfile(k(numel(k))+1:length(inputfile));
end

if plotall && (efforb>0) && nDosSites<10
    % plot all results orbital resolved, makes only sense for more than one
    % orbital
    for n=1:nDosSites
        fig(1)=figure('Position',[200, 50, 500, 300]);
        set(0,'DefaultAxesFontSize',fsz)
        plot1=plot(energy,[orbitalLDOS((n-1)*efforb+(1:efforb),:);sum(orbitalLDOS((n-1)*efforb+(1:efforb),:),1)]);
        setlabels(plot1,orb,plotrange);
        if isunix
            print_pdf(['/tmp/',inputfile,'_dx_',num2str(LDOSsites(n,1)),'_dy_',num2str(LDOSsites(n,2)),'.pdf']);
        end;
    end
end
if plotall
% impurity DOS
fig(2)=figure('Position',[200, 50, 500, 300]);
set(0,'DefaultAxesFontSize',fsz)
plot2=plot(energy,[orbitalLDOS(efforb+(1:efforb),:);sum(orbitalLDOS(efforb+(1:efforb),:),1)]);
setlabels(plot2,orb,plotrange);
if isunix
    print_pdf(['/tmp/',inputfile,'_Imp.pdf']);
end;

% NN dos
fig(3)=figure('Position',[200, 50, 500, 300]);
set(0,'DefaultAxesFontSize',fsz)
plot3=plot(energy,[orbitalLDOS(2*efforb+(1:efforb),:);sum(orbitalLDOS(2*efforb+(1:efforb),:),1)]);
setlabels(plot3,orb,plotrange);
if isunix
    print_pdf(['/tmp/',inputfile,'_Imp_NN.pdf']);
end;

% NNN dos
fig(4)=figure('Position',[200, 50, 500, 400]);
set(0,'DefaultAxesFontSize',fsz)
plot4=plot(energy,[orbitalLDOS(3*efforb+(1:efforb),:);sum(orbitalLDOS(3*efforb+(1:efforb),:),1)]);
setlabels(plot4,orb,plotrange);
if isunix
    print_pdf(['/tmp/',inputfile,'_Imp_NNN.pdf']);
end;
end;
% compare total dos
%figure5= figure('Position',[200, 50, 500, 300]);

% plot the position dependence of the summed lattice ldos
if nargin < 6
figure5= figure('Position',[150, 100, 500, 300]);
ax1 = gca;
set(0,'DefaultAxesFontSize',fsz)
else
   figure(figure5)
 hold on
end

% plot all sites or up to the plrange, or the sites given in plrange
% set up the summed ldos
if length(plrange)>1
    plrng=plrange;
else
    plrange=min(plrange,nDosSites);
    plrng=1:plrange;
end
plrange=min(plrange,nDosSites);
for n=plrng
    sumorbitalLDOS(n,:)=sum(orbitalLDOS(efforb*(n-1)+(1:efforb),:),1);
end
if isempty(plrng)
    % averaged spectra
        sumorbitalLDOS=sum(orbitalLDOS,1)*efforb/size(orbitalLDOS,1);
        % do a correction for the supercell code (when using singular_quad)
        k=strfind(inputfile,'sum');
        if isempty(k)
            % do the correction
            k1=strfind(inputfile,'_M_');
            k2=strfind(inputfile,'_ita_');
            M=str2num(inputfile(k1+3:k2-1));
        end
        sumorbitalLDOS=sumorbitalLDOS*(M^2)/(M-1)^2;
end
plot5=plot(energy,sumorbitalLDOS(plrng,:));

[e,e_minp]=find(abs(energy)==min(abs(energy)));
sumorbitalLDOS(e_minp)
% if plot_NNN
%     plot5=plot(energy,[sum(orbitalLDOS((1:efforb),:),1);sum(orbitalLDOS(efforb+(1:efforb),:),1);sum(orbitalLDOS(2*efforb+(1:efforb),:),1);sum(orbitalLDOS(3*efforb+(1:efforb),:),1)]);
% else
%     plot5=plot(energy,[sum(orbitalLDOS((1:efforb),:),1);sum(orbitalLDOS(efforb+(1:efforb),:),1);sum(orbitalLDOS(2*efforb+(1:efforb),:),1)]);
% end;

% set some special names for certain LDOSsites:
n0=0;
for n=plrng
    n0=n0+1;
    if isequal(LDOSsites(n,:),[-inf -inf])
        namestring='far away';
        lsty='-';
        lcol=[0 0 0];
    elseif isequal(LDOSsites(n,:),[0 0])
        namestring='impurity';
        lsty='--';
        lcol=coloruf2;
    elseif isequal(LDOSsites(n,:),[0 1])
        namestring='NN';
        lsty='-';
        lcol=coloruf1;
    elseif isequal(LDOSsites(n,:),[1 1])
        namestring='NNN';
        lsty='-';
        lcol=coloruf2;
    else
        namestring=['d=(',num2str(LDOSsites(n,1)),',',num2str(LDOSsites(n,2)),')'];
        lsty='';
    end
    set(plot5(n0),'DisplayName',namestring);
    if ~isempty(lsty)
       set(plot5(n0),'LineStyle',lsty,'LineWidth',2,'Color',lcol);
    end
end;
% set(plot5(1),'DisplayName','far away','LineStyle','-','LineWidth',2,'Color',[0 0 0]);
% set(plot5(2),'DisplayName','impurity','LineStyle','--','LineWidth',2,'Color',coloruf2);
% set(plot5(3),'DisplayName','NN','LineStyle','-','LineWidth',2,'Color',coloruf1);
% if plrange(end)>10
% set(plot5(4),'DisplayName','NNN','LineStyle','-','LineWidth',2,'Color',coloruf2);
% end;
xlim(plotrange);
xlabel({'\omega [eV]'});

% Create ylabel
ylabel({'LDOS [1/eV]'});

% Create legend
legend('Location','northwest')
legend('Location','northeast')
%legend show
if ~isempty(omega0)
    % put in vertical bars at the energies omega0
    %your point goes here
    l1=line([omega0 omega0],get(gca,'YLim'),'LineWidth',2,'Color',1-coloruf1);
        l2=line(-[omega0 omega0],get(gca,'YLim'),'LineWidth',2,'Color',1-coloruf2);
uistack(l1,'bottom')
uistack(l2,'bottom')
end
fndpeaks=true;
if isunix
    if isempty(plrng)
        print_pdf([inputfile_orig,'_averaged.pdf']);
    else
        print_pdf(['/tmp/',inputfile,'_tot.pdf']);
    end
    if fndpeaks
    if smoothenergy>0
        peakdistance=ceil(20*smoothenergy/de);
    else
        peakdistance=20;
    end;
    if peakdistance>200
        peakdistance=200
    end;
    energyrange=find((2*plotrange(1)<energy)+(2*plotrange(2)>energy)-1==1);
    % identify some peaks and show the positions in the plot
    try
    [imppks,implocs]=findpeaks(sum(orbitalLDOS(efforb+(1:efforb),energyrange),1),'MINPEAKDISTANCE',peakdistance,'SORTSTR','descend','NPEAKS',40);
    [NNpks,NNlocs]=findpeaks(sum(orbitalLDOS(2*efforb+(1:efforb),energyrange),1),'MINPEAKDISTANCE',peakdistance,'SORTSTR','descend','NPEAKS',40);
    [NNNpks,NNNlocs]=findpeaks(sum(orbitalLDOS(3*efforb+(1:efforb),energyrange),1),'MINPEAKDISTANCE',peakdistance,'SORTSTR','descend','NPEAKS',40);
    [farpks,farlocs]=findpeaks(sum(orbitalLDOS((1:efforb),energyrange),1),'MINPEAKDISTANCE',peakdistance,'SORTSTR','descend','NPEAKS',40);
    hold on
    energyw=energy(energyrange);
    plot(energyw(implocs),imppks,'k^','markerfacecolor',[1 0 0])
    plot(energyw(NNlocs),NNpks,'kv','markerfacecolor',coloruf1)
    plot(energyw(NNNlocs),NNNpks,'k*','markerfacecolor',coloruf2)
    plot(energyw(farlocs),farpks,'ko','markerfacecolor',[0 1 0])
    xlim(2*plotrange);
    f={[energyw(implocs);imppks],[energyw(NNlocs);NNpks],[energyw(NNNlocs);NNNpks],[energyw(farlocs);farpks]};
    catch
        disp('problem with findpeaks, probably no Signal Processing Toolbox available');
    end
    % give back the peak positions of the 5 largest peaks
    end;
end;

    if fndpeaks

%peak detection plot 
figure6= figure('Position',[200, 50, 500, 300]);
data=[sum(orbitalLDOSImpNNN)./sum(orbitalLDOS(2*efforb+(1:efforb),:),1);1./sum(orbitalLDOS(3*efforb+(1:efforb),:),1).*sum(orbitalLDOS(2*efforb+(1:efforb),:),1)];
plot6=plot(energy,data);
set(plot6(1),'DisplayName','NNN/NN','Color',coloruf2);
set(plot6(2),'DisplayName','NN/NNN','Color',coloruf1);
xlim(plotrange);
ylim([0,max(data(:))]);
xlabel({'\omega'});

% Create ylabel
ylabel({'LDOS rel'});

% Create legend
legend1 = legend(gca,'show');
set(legend1,'Location','Best');
if isunix
  %  print_pdf(['/tmp/',inputfile,'_rel.pdf']);
end;
    end
end

function setlabels(plot,orb,range)
szorb=size(orb,2);
for i=1:szorb
set(plot(i),'DisplayName',orb{i});
end;
xlim(range);
% Create xlabel
xlabel({'\omega [eV]'});

% Create ylabel
ylabel({'LDOS [1/eV]'});

% Create legend
legend show
end
