function h=plot_homogeneous_dos_v2(inputfile,smoothenergy)

% Modified homogeneous_dos.m
% takes \Delta_ij as input and constructs \Delta_i0.

if (~exist('plotrange','var'))
    plotrange=[-0.2 0.2];
end;

if (~exist('tetra','var'))
    tetra=false;
end;

if nargin <1
    % load relevant files
    TB_file='TB_hamiltonian_FeSe_2D.mat'
    %latticeVectors = latticeVector;
    BdGfileName='BdG_homogeneous_FeSe_Milan_GammaCut_2_N_9(1).mat'
    Gamma_file='Gamma_FeSe_Milan_GammaCut_2.mat'
    M = input('Enter no of k-points   ');% no of K points in x
    ita = input('Enter ita   '); % broadening
    firstEnergy = input('Enter starting energy   ');
    lastEnergy = input('Enter last energy   ');
    nEnergyPoints = input('Enter no of energy points   ');
else
    % load the input file to set the variables, gave up the old .mat file
    % format
     read_input_file=inputfile;
     read_input;
     read_input_file
end;

if nargin < 2
    if ~tetra
        smoothenergy=0;
    else
        %default: smoothing of 2 meV if "tetra" is set
        smoothenergy=0.002;
    end
end
sqstring='';
if (exist('singular_quad','var'))
    sqstring='sum';
end;
if ~exist('sublattice','var')
    sublattice=1;
end;
% special convention for homogeneous DOS
M=N*M;
LDOSfileName0 = [casestring,'_Vimp_', num2str(Vimp),  '_N_', num2str(N)];

if ~tetra
    LDOSfileName = [LDOSfileName0,sqstring , '_M_', num2str(M),'_ita_', num2str(ita)];
else
    LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_tetra_corr1']
end;

disp('reading k-space calculated DOS ...');
load(LDOSfileName,'-mat')

% some smoothing if necessary
de=energy(2)-energy(1);
smooth=floor(smoothenergy/de);
bandDOSNormal=sg_smooth(bandDOSNormal,smooth);
bandDOS=sg_smooth(bandDOS,smooth);


totalDOSNormal = sum(bandDOSNormal,1);
totalDOS = sum(bandDOS,1);

k = findstr(inputfile, '/');
if ~isempty(k)
    	inputfile=inputfile(k(numel(k))+1:length(inputfile));
end;

fsz=14;

% Plotting
fig1= figure('Position',[150, 100, 500, 300]);
set(0,'DefaultAxesFontSize',fsz)
%hold on
if sublattice==0
    sublatticefactor=1;
else
    sublatticefactor=0.5;
end;
plotres=[sublatticefactor*totalDOSNormal;sublatticefactor*totalDOS];
plot1=plot(energy,plotres);
set(plot1(1),'LineWidth',2,'LineStyle','-','DisplayName','Normal state','Color',[0 0 0]);
set(plot1(2),'Color',[1 0 0],'DisplayName','SC state');
xlim(plotrange);
ylim_curr = get(gca,'ylim');
set(gca, 'ylim', [0 ylim_curr(2)]);
%,'k');
%hold;
%plot(energy,(5/nOrbitals)*totalDOS, 'r');
% Create xlabel
xlabel({'\omega [eV]'});

% Create ylabel
ylabel({'DOS [1/eV]'});
%axis('square'); title('Normal Vs SC dos')
% Create legend
legend show

if tetra
    pdffile1=['/tmp/',inputfile,'_normal_SC_tetra',num2str(smoothenergy),'.pdf']
else
    pdffile1=['/tmp/',inputfile,'_normal_SC.pdf']
end

print_pdf(pdffile1);

coloruf1=[250 	70 	22 ]/255;
coloruf2= [0 	48 	135]/255;

fig2= figure('Position',[200, 50, 500, 300]);
set(0,'DefaultAxesFontSize',fsz)
totDOS=sublatticefactor*totalDOS;
%plot(energy, (5/nOrbitals)*totalDOS, 'k');
plot1=plot(energy, [totDOS;sublatticefactor*totalDOSNormal;bandDOS(1:sublatticefactor*nOrbitals,:)]);
%, 'r');
% plot(energy, bandDOS(2,:), 'g');
% plot(energy, bandDOS(3,:), 'c');
% plot(energy, bandDOS(4,:), 'm');
% plot(energy, bandDOS(5,:), 'b');

set(plot1(1),'LineWidth',2,'DisplayName','total',...
    'Color',[0 0 0]);
set(plot1(2),'LineWidth',2,'LineStyle',':','Color',[0 0 0],'DisplayName','normal state');
linestyles={':','-.','-','--','-'};
colors={[1 0 0], [1 0 0], [0 1 0], [1 0 0], [0 0 1], [0 0 0]};
markers={'','','o','','',''};
Displaynames={'d_{z^2}','d_{x^2-y^2}','d_{yz}','d_{xz}','d_{xy}'};
for n=1:sublatticefactor*nOrbitals
    if isempty(markers{n})
        set(plot1(n+2),'LineStyle',linestyles{n},'Color',colors{n},'DisplayName',Displaynames{n});
    else
        set(plot1(n+2),'LineStyle',linestyles{n},'Color',colors{n},'DisplayName',Displaynames{n},'Marker',markers{n});
    end
end;
%set(plot1(3),'LineStyle',':','Color',[1 0 0],'DisplayName','d_{z^2}');
%set(plot1(4),'LineStyle','-.','Color',[1 0 0],'DisplayName','d_{x^2-y^2}');
%set(plot1(5),'MarkerSize',2,'Marker','o','Color',[0 1 0],...
%    'DisplayName','d_{yz}');
%set(plot1(6),'LineStyle','--','Color',[1 0 0],'DisplayName','d_{xz}');
%set(plot1(7),'Color',[0 0 1],'DisplayName','d_{xy}');
%set(plot1(2),'LineWidth',2,'LineStyle',':','Color',[0 0 0],'DisplayName','normal state');

xlim(plotrange);
ylim_curr = get(gca,'ylim');
set(gca, 'ylim', [0 ylim_curr(2)]);
%axis('square'); 
%title('Orbital resolved SC dos')
% Create xlabel
xlabel({'\omega [eV]'});

% Create ylabel
ylabel({'DOS [1/eV]'});
% Create legend
legend show
if tetra
    pdffile2=['/tmp/',inputfile,'_orbital_SC_tetra',num2str(smoothenergy),'.pdf']
else
    pdffile2=['/tmp/',inputfile,sqstring,'_orbital_SC.pdf']
end
print_pdf(pdffile2);
