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
    if isempty(strfind(inputfile, 'tetra'))
        smoothenergy=0;
    else
        %default: smoothing of 2 meV if "tetra" in filename
        smoothenergy=0.002;
    end
end


% special convention for homogeneous DOS
M=N*M;
LDOSfileName0 = [casestring,'_Vimp_', num2str(Vimp),  '_N_', num2str(N)];

if ~tetra
    LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_ita_', num2str(ita)];
else
    LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_tetra_corr']
end;

disp('reading k-space calculated DOS ...');
load(LDOSfileName,'-mat')

% some smoothing if necessary
de=energy(2)-energy(1);
smooth=floor(smoothenergy/de);
bandDOSNormal=sg_smooth(bandDOSNormal,smooth);
bandDOS=sg_smooth(bandDOS,smooth);


totalDOSNormal = sum(bandDOSNormal);
totalDOS = sum(bandDOS);

k = findstr(inputfile, '/');
if ~isempty(k)
    	inputfile=inputfile(k(numel(k))+1:length(inputfile));
end;


% Plotting
fig1=figure; 
plot1=plot(energy,[(5/nOrbitals)*totalDOSNormal;(5/nOrbitals)*totalDOS]);
set(plot1(1),'DisplayName','Normal state','Color',[0 0 0]);
set(plot1(2),'Color',[1 0 0],'DisplayName','SC state');
%,'k');
%hold;
%plot(energy,(5/nOrbitals)*totalDOS, 'r');
% Create xlabel
xlabel({'\omega'});

% Create ylabel
ylabel({'DOS [1/eV]'});
axis('square'); title('Normal Vs SC dos')
% Create legend
legend show
xlim(plotrange);

print_pdf(['/tmp/',inputfile,'_normal_SC.pdf']);

fig2=figure;
totDOS=(5/nOrbitals)*totalDOS;
%plot(energy, (5/nOrbitals)*totalDOS, 'k');
plot1=plot(energy, [totDOS;bandDOS(1:5,:)]);
%, 'r');
% plot(energy, bandDOS(2,:), 'g');
% plot(energy, bandDOS(3,:), 'c');
% plot(energy, bandDOS(4,:), 'm');
% plot(energy, bandDOS(5,:), 'b');

set(plot1(1),'DisplayName','total','Color',[0 0 0]);
set(plot1(2),'Color',[1 0 0],'DisplayName','orbital1');
set(plot1(3),'Color',[0 1 0],'DisplayName','orbital2');
set(plot1(4),'Color',[0 1 1],'DisplayName','orbital3');
set(plot1(5),'Color',[1 0 1],'DisplayName','orbital4');
set(plot1(6),'Color',[0 0 1],'DisplayName','orbital5');

axis('square'); title('Orbital resolved SC dos')
% Create xlabel
xlabel({'\omega [eV]'});

% Create ylabel
ylabel({'DOS [1/eV]'});
% Create legend
legend show
xlim(plotrange);
print_pdf(['/tmp/',inputfile,'_orbital_DOS_SC.pdf']);
