function h=plot_k_space_gap(inputfile)
    % load the input file to set the variables, gave up the old .mat file
    % format
     read_input_file=inputfile;
     read_input;
     read_input_file
if (~exist('Vimp','var'))
    if abs(Vimp)>0
        disp('Warning: finite impurity potential, not homogeneous case.')
    end;
end;
load(TB_file,'-mat');
nOrbitals = size(TBparameters,1);
    load(Gamma_file,'-mat');
    load(BdGfileName,'-mat');
N = sqrt(size(delta,1)/nOrbitals);
nUnitCellsDelta = size(latticeVectorsSC,1);
deltaCenter = zeros(nOrbitals, nOrbitals, size(latticeVectorsSC,1));
jCell = [ceil(N/2) ceil(N/2)];
for i = 1:nUnitCellsDelta
    iCell = jCell + latticeVectorsSC(i,:);
    [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, iCell, jCell);
    deltaCenter(:,:,i) = delta(iRange, jRange);
end
delta = deltaCenter;



nUnitCells = size(latticeVector,1);
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) = ...
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) - mu*eye(nOrbitals);

M=N*M;
M=100;
kx = (2*pi/M)*(0:(M - 1)) + pi/M;
ky = kx;
delKx = kx(2)-kx(1);
delKy = delKx;
%energy = linspace(firstEnergy, lastEnergy, nEnergyPoints);
kSpaceEigenValues = zeros(M, M, 2*nOrbitals);
kSpaceEigenVectors = zeros(M, M, 2*nOrbitals, 2*nOrbitals);
kSpaceEigenValuesNormal = zeros(M, M, nOrbitals);
kSpacegapall = zeros(M, M, nOrbitals,nOrbitals);
kSpaceEigenVectorsNormal = zeros(M, M, nOrbitals, nOrbitals);

% some double code from homogeneous_dos_v2.m !
for iKx = 1:M
        for iKy = 1:M
            k = [kx(iKx) ky(iKy)];
            % diagonalizing for normal state DOS
            kSpaceHopping = zeros(nOrbitals,nOrbitals);
            for iUnitCell = 1:nUnitCells
                iLatticeVector = latticeVector(iUnitCell,:);
                kSpaceHopping = kSpaceHopping + TBparameters(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
            end                        
            [eigVectorNormal eigValueNormal] = eig(kSpaceHopping);
            [eigValueKNormal sortingIndexNormal] = sort(real(diag(eigValueNormal)));
            eigVectorKNormal = (eigVectorNormal(:,sortingIndexNormal))';
            kSpaceEigenValuesNormal(iKx, iKy, :) = eigValueKNormal;
            kSpaceEigenVectorsNormal(iKx, iKy, :,:) =  eigVectorKNormal;
            % diagonalizing for SC state DOS
            kSpaceGap =zeros(nOrbitals,nOrbitals);
            for iUnitCellDelta = 1:nUnitCellsDelta
                iLatticeVectorDelta = latticeVectorsSC(iUnitCellDelta,:);
                kSpaceGap = kSpaceGap + delta(:,:,iUnitCellDelta)*exp(1i*(iLatticeVectorDelta*k'));
            end
            kSpacegapall(iKx,iKy,:,:)=kSpaceGap(:,:);
            kSpaceHamiltonian = [kSpaceHopping -kSpaceGap; -kSpaceGap' -kSpaceHopping];
            [eigVector eigValue] = eig(kSpaceHamiltonian);
            [eigValueK sortingIndex] = sort(real(diag(eigValue)));
            eigVectorK = (eigVector(:,sortingIndex))';
            kSpaceEigenValues(iKx, iKy, :) = eigValueK;
            kSpaceEigenVectors(iKx, iKy, :,:) =  eigVectorK;
        end
        if  mod(iKx,10)==0   
            disp(['Done ',num2str(iKx), ' of ',num2str(M),' kx values.']);
        end;
end
 surf(real(kSpacegapall(:,:,1,1)))

