function f=make_plot_impurity_dos(inputfile,smoothenergy)
if nargin <1
    inputfile='LDOS_FeSe_Milan_Gamma_Vimp_4_N_9_M_40_ita_0.003'
end;
plotall=false;
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

% load input
load(inputfile,'-mat')
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
plot1=plot(energy,[orbitalLDOSFarAway;sum(orbitalLDOSFarAway)]);
setlabels(plot1,orb,plotrange);

k = findstr(inputfile, '/');
if ~isempty(k)
    	inputfile=inputfile(k(numel(k))+1:length(inputfile));
end;

print_pdf(['/tmp/',inputfile,'_far_away.pdf']);
if plotall
% impurity DOS
figure2=figure('Position',[200, 50, 500, 300]);
plot2=plot(energy,[orbitalLDOSImp;sum(orbitalLDOSImp)]);
setlabels(plot2,orb,plotrange);
print_pdf(['/tmp/',inputfile,'_Imp.pdf']);

% NN dos
figure3=figure('Position',[200, 50, 500, 300]);
plot3=plot(energy,[orbitalLDOSImpNN;sum(orbitalLDOSImpNN)]);
setlabels(plot3,orb,plotrange);
print_pdf(['/tmp/',inputfile,'_Imp_NN.pdf']);

% NNN dos
figure4=figure('Position',[200, 50, 500, 300]);
plot4=plot(energy,[orbitalLDOSImpNNN;sum(orbitalLDOSImpNNN)]);
setlabels(plot4,orb,plotrange);
print_pdf(['/tmp/',inputfile,'_Imp_NNN.pdf']);
end;
% compare total dos
figure5= figure('Position',[200, 50, 500, 300]);

plot5=plot(energy,[sum(orbitalLDOSFarAway);sum(orbitalLDOSImp);sum(orbitalLDOSImpNN);sum(orbitalLDOSImpNNN)]);
set(plot5(1),'DisplayName','tot far away','LineStyle','-','LineWidth',2,'Color',[0 0 0]);
set(plot5(2),'DisplayName','tot impurity','LineStyle','--','LineWidth',2,'Color',[0 0 0]);
set(plot5(3),'DisplayName','tot NN','LineStyle','-','LineWidth',2,'Color',coloruf1);
set(plot5(4),'DisplayName','tot NNN','LineStyle','-','LineWidth',2,'Color',coloruf2);
xlim(plotrange);
xlabel({'\omega'});

% Create ylabel
ylabel({'DOS [1/eV]'});

% Create legend
legend show
print_pdf(['/tmp/',inputfile,'_tot.pdf']);

end

function setlabels(plot,orb,range)
for i=1:6
set(plot(i),'DisplayName',orb{i});
end;
xlim(range);
% Create xlabel
xlabel({'\omega [eV]'});

% Create ylabel
ylabel({'DOS [1/eV]'});

% Create legend
legend show
end
