function [kSpaceEigenValuesNormal,kSpaceEigenVectorsNormal]=bands_klist(inputfile,klistfile)

% Modified homogeneous_dos.m
% takes \Delta_ij as input and constructs \Delta_i0.

if nargin <1
    % load relevant files
    TB_file='TB_hamiltonian_FeSe_2D.mat'
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

load(TB_file,'-mat');
nOrbitals = size(TBparameters,1);
    mu=0;


nUnitCells = size(latticeVector,1);
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) = ...
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) - mu*eye(nOrbitals);
% fix for the lattice vectors
latticeVector=latticeVector(:,1:2);

% read in k values
kin=load(klistfile);
szk=size(kin);
kSpaceEigenValuesNormal = zeros(szk(1), nOrbitals);
kSpaceEigenVectorsNormal = zeros(szk(1), nOrbitals, nOrbitals);
for iKx = 1:szk(1)
            k = [kin(iKx,1) kin(iKx,2)];
            % diagonalizing for normal state DOS
            kSpaceHopping = zeros(nOrbitals,nOrbitals);
            for iUnitCell = 1:nUnitCells
                iLatticeVector = latticeVector(iUnitCell,:);
                kSpaceHopping = kSpaceHopping + TBparameters(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
            end                        
            [eigVectorNormal, eigValueNormal] = eig(kSpaceHopping);
            [eigValueKNormal, sortingIndexNormal] = sort(real(diag(eigValueNormal)));
            eigVectorKNormal = (eigVectorNormal(:,sortingIndexNormal))';
            kSpaceEigenValuesNormal(iKx, :) = eigValueKNormal;
            kSpaceEigenVectorsNormal(iKx, :,:) =  eigVectorKNormal;
end


