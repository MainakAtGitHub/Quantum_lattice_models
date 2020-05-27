function [kSpaceEigenValuesNormal,kSpaceEigenVectorsNormal]=bands_precalc_klist(inputfile,klist)

% for given inputfile (*.txt) and an array containing the momentum space Hamiltonian, calculate the eigenvalues and eigenvectors
% of the band structure

% call the function that does the Fourier transform
[kSpaceHamiltonian]=Hk_klist(inputfile,klist);

% by default set the chemical potential to zero
 mu=1.062;
    % load the input file to set the variables, gave up the old .mat file
    % format
     read_input_file=inputfile;
     read_input;
     read_input_file


% load the hoppings from the file
load(TB_file,'-mat');
nOrbitals = size(TBparameters,1);

% number of hoppings in band structure
nUnitCells = size(latticeVector,1);

% allow to set the chemical potential
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) = ...
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) - mu*eye(nOrbitals);

% fix for the lattice vectors (no 3rd component allowed)
latticeVector=latticeVector(:,1:2);

% check the size of the input array
szk=size(kSpaceHamiltonian);
if ~(szk(2)==nOrbitals)
    disp('mismatch of number of orbitals, check input')
end
% initialize the result variables
kSpaceEigenValuesNormal = zeros(szk(1), nOrbitals);
kSpaceEigenVectorsNormal = zeros(szk(1), nOrbitals, nOrbitals);

%iterate over the momentum vectors
for iKx = 1:szk(1)
	    % diagonalize the Hamiltonian
            [eigVectorNormal, eigValueNormal] = eig(squeeze(kSpaceHamiltonian(iKx,:,:)));
            % sort the eigenvalues and give back the order
            [eigValueKNormal, sortingIndexNormal] = sort(real(diag(eigValueNormal)));
	    % put eigenvectors in the corresponding order
            eigVectorKNormal = (eigVectorNormal(:,sortingIndexNormal))';
	    % put them into the restult array
            kSpaceEigenValuesNormal(iKx, :) = eigValueKNormal;
            kSpaceEigenVectorsNormal(iKx, :,:) =  eigVectorKNormal;
end
kSpaceEigenValuesNormal=kSpaceEigenValuesNormal-mu;