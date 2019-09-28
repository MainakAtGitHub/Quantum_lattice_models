function [kspace_gap]=Deltak_klist(inputfile,klist)

% read in from a (converged) calculation, the real space gaps and do the
% Fourier transformation to momentum space, giving back an array with the
% order parameter in momentum space for each of the k-values

    % load the input file to set the variables, gave up the old .mat file
    % format
     read_input_file=inputfile;
     % call the script (not function that parses the input file and sets up
     % the correct variables)
     read_input;
    % read_input_file



% load variables from the Gamma file and the (converged) BdGfile:
load(Gamma_file,'-mat');
    load(BdGfileName,'-mat');
    % get the number of sites (does not work for non-square systems yet,
    % please be careful)
N = sqrt(size(delta,1)/nOrbitals);
% number of NN vectors from the pairing interaction
nUnitCellsDelta = size(latticeVectorsSC,1);
% set up array to hold the real space pairing terms
deltaCenter = zeros(nOrbitals, nOrbitals, size(latticeVectorsSC,1));
% Message if the data is from a calculation with impurity
if ~(Vimp==0)
    disp('Seems to be a calculation with impurity, results are not reliable!')
end
% arbitrarily select a center point (should be independend for a converged
% homogeneous system)
jCell = [ceil(N/2) ceil(N/2)];
for i = 1:nUnitCellsDelta
    iCell = jCell + latticeVectorsSC(i,:);
    [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, iCell, jCell);
    deltaCenter(:,:,i) = delta(iRange, jRange);
end
delta = deltaCenter;



% read in k values if klist is a string, else just use the klist as data
if isa(klist,'char') || isa(klist,'string') 
	kin=load(klist);
else
    kin=klist;
end
szk=size(kin);
% set up result array
kspace_gap = zeros(szk(1), nOrbitals, nOrbitals);
% iterate over k-values
for iK = 1:szk
            k = [kin(iK,1) kin(iK,2)];
            % initialize zero matrix
            kSpaceGap = zeros(nOrbitals, nOrbitals);
            % do the Fourier transform
            for iUnitCellDelta = 1:nUnitCellsDelta
                iLatticeVectorDelta = latticeVectorsSC(iUnitCellDelta,:);
                kSpaceGap = kSpaceGap + delta(:,:,iUnitCellDelta)*exp(1i*(iLatticeVectorDelta*k'));
            end
            % put in the result into array
            kspace_gap(iK, :,:)=kSpaceGap;
end
