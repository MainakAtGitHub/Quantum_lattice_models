function r=write_input_impurity_dos(outputfile)

% calculate the impurity density of states
% inputfile : file that contains the parameters
% division : divide task into division parts
% part : calculate this part part = 1...division

% if no inputfile is specified, read in from the command line
if nargin <1
    outputfile= 'standart_input_imp_dos.mat';
end;

    % input parameters
    N = 13;
    Vimp = 0.4;
    M = 20;% no of K points in x
    ita = 0.005; % broadening
    firstEnergy = -1;
    lastEnergy = 1;
    nEnergyPoints = 2500;
    % set some input filenames
    % tight binding model
    TB_file='TB_hamiltonian_FeSe_2D.mat';
    % interactions
    % old input
    % BdGfileName = ['BdG_Impurity_FeSe', '_N_', num2str(N),'_Vimp_', num2str(Vimp)];
    % casestring='LDOS_FeSe_Milan_Gamma';
    % load(TB_file);  % not necessary to save again
    % possibly not necessary?
    % latticeVectors = latticeVector;
    % save variables for next run
    
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
end;
    save(outputfile)
