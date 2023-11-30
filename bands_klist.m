function [kSpaceEigenValuesNormal,kSpaceEigenVectorsNormal]=bands_klist(inputfile,klist)

% for given inputfile (*.txt) and a file containing k-points, calculate the eigenvalues and eigenvectors
% of the band structure

% by default set the chemical potential to zero
 mu=0;
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

% read in k values if klist is a string, else just use the klist as data
if isa(klist,'char') || isa(klist,'string')
	kin=load(klist);
else
    kin=klist;
end
szk=size(kin);
% initialize the result variables
kSpaceEigenValuesNormal = zeros(szk(1), nOrbitals);
kSpaceEigenVectorsNormal = zeros(szk(1), nOrbitals, nOrbitals);

%iterate over the momentum vectors
for iKx = 1:szk(1)
	    % set the momentum vector (in 2D)
            k = [kin(iKx,1) kin(iKx,2)];
	    % initialize the Hamiltonian to zeros
            kSpaceHopping = zeros(nOrbitals,nOrbitals);
	    % add all terms coming from each hopping
            for iUnitCell = 1:nUnitCells
                iLatticeVector = latticeVector(iUnitCell,:);
		% do (manually) the Fourier transform and add to the Hamiltonian
                kSpaceHopping = kSpaceHopping + TBparameters(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
            end
	    % diagonalize the Hamiltonian
            [eigVectorNormal, eigValueNormal] = eig(kSpaceHopping);
            [eigValueKNormal, sortingIndexNormal] = sort(real(diag(eigValueNormal)));
	    % sort the eigenvalues and eigenvectors
            eigVectorKNormal = (eigVectorNormal(:,sortingIndexNormal))';
	    % put them into the restult array
            kSpaceEigenValuesNormal(iKx, :) = eigValueKNormal;
            kSpaceEigenVectorsNormal(iKx, :,:) =  eigVectorKNormal;
end
