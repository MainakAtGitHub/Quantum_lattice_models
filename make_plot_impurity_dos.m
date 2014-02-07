function f=make_plot_impurity_dos(inputfile,smoothenergy,plotall,omega0)
if nargin <1
    inputfile='LDOS_FeSe_Milan_Gamma_Vimp_4_N_9_M_40_ita_0.003'
end;
if nargin <3
    plotall=false;
end;
coloruf1=[250 	70 	22 ]/255;
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

if (~exist('plotrange','var'))
    plotrange=[-0.2 0.2];
end;

orb={'orbital1','orbital2','orbital3','orbital4','orbital5','total'};
fsz=14;
% load input
try
    load(inputfile,'-mat')
catch
    load(inputfile)
end;
% some smoothing if necessary
de=energy(2)-energy(1);
smooth=floor(smoothenergy/de);
orbitalLDOSFarAway=sg_smooth(orbitalLDOSFarAway,smooth);
orbitalLDOSImp=sg_smooth(orbitalLDOSImp,smooth);
orbitalLDOSImpNN=sg_smooth(orbitalLDOSImpNN,smooth);
orbitalLDOSImpNNN=sg_smooth(orbitalLDOSImpNNN,smooth);
% fix negative values (numerical error, smoothing artefacts)
orbitalLDOSFarAway(orbitalLDOSFarAway<0)=0;
orbitalLDOSImp(orbitalLDOSImp<0)=0;
orbitalLDOSImpNN(orbitalLDOSImpNN<0)=0;
orbitalLDOSImpNNN(orbitalLDOSImpNNN<0)=0;
% DOS far away (without impurity)
figure1=figure('Position',[200, 50, 500, 300]);
set(0,'DefaultAxesFontSize',fsz)
plot1=plot(energy,[orbitalLDOSFarAway;sum(orbitalLDOSFarAway)]);
setlabels(plot1,orb,plotrange);

k = findstr(inputfile, '/');
if ~isempty(k)
    	inputfile=inputfile(k(numel(k))+1:length(inputfile));
end;

if plotall
if isunix
    print_pdf(['/tmp/',inputfile,'_far_away.pdf']);
end;
% impurity DOS
figure2=figure('Position',[200, 50, 500, 300]);
set(0,'DefaultAxesFontSize',fsz)
plot2=plot(energy,[orbitalLDOSImp;sum(orbitalLDOSImp)]);
setlabels(plot2,orb,plotrange);
if isunix
    print_pdf(['/tmp/',inputfile,'_Imp.pdf']);
end;

% NN dos
figure3=figure('Position',[200, 50, 500, 300]);
set(0,'DefaultAxesFontSize',fsz)
plot3=plot(energy,[orbitalLDOSImpNN;sum(orbitalLDOSImpNN)]);
setlabels(plot3,orb,plotrange);
if isunix
    print_pdf(['/tmp/',inputfile,'_Imp_NN.pdf']);
end;

% NNN dos
figure4=figure('Position',[200, 50, 500, 300]);
set(0,'DefaultAxesFontSize',fsz)
plot4=plot(energy,[orbitalLDOSImpNNN;sum(orbitalLDOSImpNNN)]);
setlabels(plot4,orb,plotrange);
if isunix
    print_pdf(['/tmp/',inputfile,'_Imp_NNN.pdf']);
end;
end;
% compare total dos
%figure5= figure('Position',[200, 50, 500, 300]);
figure5= figure('Position',[150, 100, 500, 300]);
set(0,'DefaultAxesFontSize',fsz)

plot5=plot(energy,[sum(orbitalLDOSFarAway);sum(orbitalLDOSImp);sum(orbitalLDOSImpNN);sum(orbitalLDOSImpNNN)]);
set(plot5(1),'DisplayName','tot far away','LineStyle','-','LineWidth',2,'Color',[0 0 0]);
set(plot5(2),'DisplayName','tot impurity','LineStyle','--','LineWidth',2,'Color',[0 0 0]);
set(plot5(3),'DisplayName','tot NN','LineStyle','-','LineWidth',2,'Color',coloruf1);
set(plot5(4),'DisplayName','tot NNN','LineStyle','-','LineWidth',2,'Color',coloruf2);
xlim(plotrange);
xlabel({'\omega [eV]'});

% Create ylabel
ylabel({'LDOS [1/eV]'});

% Create legend
legend show
if exist('omega0','var')
    % put in vertical bars at the energies omega0
    %your point goes here
    l1=line([omega0 omega0],get(gca,'YLim'),'LineWidth',2,'Color',1-coloruf1);
        l2=line(-[omega0 omega0],get(gca,'YLim'),'LineWidth',2,'Color',1-coloruf2);
uistack(l1,'bottom')
uistack(l2,'bottom')
end;
if isunix
    print_pdf(['/tmp/',inputfile,'_tot.pdf']);
    if smoothenergy>0
        peakdistance=ceil(20*smoothenergy/de);
    else
        peakdistance=20;
    end;
    energyrange=find((2*plotrange(1)<energy)+(2*plotrange(2)>energy)-1==1);
    % identify some peaks and show the positions in the plot
    [imppks,implocs]=findpeaks(sum(orbitalLDOSImp(:,energyrange)),'MINPEAKDISTANCE',peakdistance,'SORTSTR','descend','NPEAKS',40);
    [NNpks,NNlocs]=findpeaks(sum(orbitalLDOSImpNN(:,energyrange)),'MINPEAKDISTANCE',peakdistance,'SORTSTR','descend','NPEAKS',40);
    [NNNpks,NNNlocs]=findpeaks(sum(orbitalLDOSImpNNN(:,energyrange)),'MINPEAKDISTANCE',peakdistance,'SORTSTR','descend','NPEAKS',40);
    [farpks,farlocs]=findpeaks(sum(orbitalLDOSFarAway(:,energyrange)),'MINPEAKDISTANCE',peakdistance,'SORTSTR','descend','NPEAKS',40);
    hold on
    energyw=energy(energyrange);
    plot(energyw(implocs),imppks,'k^','markerfacecolor',[1 0 0])
    plot(energyw(NNlocs),NNpks,'kv','markerfacecolor',coloruf1)
    plot(energyw(NNNlocs),NNNpks,'k*','markerfacecolor',coloruf2)
    plot(energyw(farlocs),farpks,'ko','markerfacecolor',[0 1 0])
    xlim(2*plotrange);
    f={[energyw(implocs);imppks],[energyw(NNlocs);NNpks],[energyw(NNNlocs);NNNpks],[energyw(farlocs);farpks]};
    % give back the peak positions of the 5 largest peaks
end;


%peak detection plot 
figure6= figure('Position',[200, 50, 500, 300]);
data=[sum(orbitalLDOSImpNNN)./sum(orbitalLDOSImpNN);1./sum(orbitalLDOSImpNNN).*sum(orbitalLDOSImpNN)];
plot6=plot(energy,data);
set(plot6(1),'DisplayName','NNN/NN','Color',coloruf2);
set(plot6(2),'DisplayName','NN/NNN','Color',coloruf1);
xlim(2*plotrange);
ylim([0,max(data(:))]);
xlabel({'\omega'});

% Create ylabel
ylabel({'LDOS rel'});

% Create legend
legend1 = legend(gca,'show');
set(legend1,'Location','Best');
if isunix
    print_pdf(['/tmp/',inputfile,'_rel.pdf']);
end;

end

function setlabels(plot,orb,range)
for i=1:6
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
