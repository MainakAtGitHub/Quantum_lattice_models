function r=write_input_impurity_dos(inputfile)

% write the input file for the calculations
% outputfile : file to be written (can be omitted)

% if no inputfile is specified, read in from the command line
if nargin <1
    outputfile= 'standart_input_imp_dos.mat';
end;

% input parameters (modify as needed!)
N = 9
Vimp=0
%Vimp = 0.4
%Vimp = 0.8;
%Vimp = 1.6;
%Vimp = 3.2;
M = 20;% no of K points in x
ita = 0.001; % broadening
firstEnergy = -1
lastEnergy = 1
nEnergyPoints = 2500
% set some input filenames
% tight binding model
TB_file='TB_hamiltonian_FeSe_2D.mat'
% interactions
% old input
% BdGfileName = ['BdG_Impurity_FeSe', '_N_', num2str(N),'_Vimp_', num2str(Vimp)];
% casestring='LDOS_FeSe_Milan_Gamma';
% load(TB_file);  % not necessary to save again
% possibly not necessary?
% latticeVectors = latticeVector;
% save variables for next run
alpha = .25;%input('enter alpha    '); % self-consistency parametser
beta1 = .7;%input('enter beta1    ');
beta2 = .9;%input('enter beta1    ');
% suggestion by Maria take <beta=0.85> (take larger values for larger
% Vimp!)
%beta1 = .88;%input('enter beta1    ');
%beta2 = .95;%input('enter beta1    ');
deltaTol = 1e-4;%input('enter tolerance for gap convergence    ');
deltaTol = 1e-6;%input('enter tolerance for gap convergence    ');
nTol = 1e-3;%input('enter tolerance for electron density convergence    ');
nTol = 1e-5;%input('enter tolerance for electron density convergence    ');
maxLoop = 500;%input('enter maxloop     '); % max no of iterations for self consistency
nOrbitals = 10;
n0 = 1.2*nOrbitals; % no. of valence electrons per unit cell
kT = .01;
casename='milan';
casename='tom';
switch casename
    case 'milan'
        Gamma_file='Gamma_FeSe_10_orbital_Milan_symmetrized.mat';
        casestring='LDOS_FeSe_Milan_';
        BdGfileName = 'BdG_Impurity_FeSe_N_9_Vimp_4.mat';
    case 'tom'
        casestring='LDOS_FeSe_Tom_';
        BdGfileName = ['BdG_Impurity_FeSe_Toms_BS_6Dec13_N_' num2str(N),'_Vimp_', num2str(Vimp)];
        BdGfileName = 'BdG_Impurity_FeSe_Toms_BS_6Dec13_N_13_Vimp_0.4.mat'
        %  BdGfileName ='BdG_Impurity_FeSe_Toms_BS_6Dec13_N_13_Vimp_0_no_SC.mat';
        Gamma_file='Gamma_FeSe_Toms_BS_6Dec13_cut_2.mat'
        Gamma_file='Gamma_FeSe_Toms_BS_6Dec13_GammaCut_3.mat'
        BdGfileName = ['BdG_Impurity_FeSe_Toms_BS_6Dec13',num2str(maxLoop), '_N_', num2str(N),'_Vimp_', num2str(Vimp)];
        BdGfileName = ['BdG_Imp_FeSe_new',num2str(maxLoop), '_N_', num2str(N),'_Vimp_', num2str(Vimp)]
        BdGfileName = ['initial_guess_N_9_l1.mat']
        BdGfileName = ['initial_guess_N_9_complex.mat']
        input_fileName = BdGfileName;%['BdG_homogeneous_FeSe_Toms_BS_6Dec13', '_N_', num2str(N),'_GammaCut_',num2str(2),'.mat'];
end;
save(inputfile)
