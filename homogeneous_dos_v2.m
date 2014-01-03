% Modified homogeneous_dos.m
% takes \Delta_ij as input and constructs \Delta_i0.


% load relevant files
load TB_hamiltonian_FeSe_2D.mat
latticeVectors = latticeVector;
load BdG_homogeneous_FeSe_Milan_GammaCut_2_N_9(1).mat
load Gamma_FeSe_Milan_GammaCut_2.mat
M = input('Enter no of k-points   ');% no of K points in x
ita = input('Enter ita   '); % broadening
firstEnergy = input('Enter starting energy   ');
lastEnergy = input('Enter last energy   ');
nEnergyPoints = input('Enter no of enery points   ');


nOrbitals = size(TBparameters,1);
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


nUnitCells = size(latticeVectors,1);
TBparameters(:,:,(latticeVectors(:,1)==0) & (latticeVectors(:,2)==0)) = ...
TBparameters(:,:,(latticeVectors(:,1)==0) & (latticeVectors(:,2)==0)) - mu*eye(nOrbitals);


kx = (2*pi/M)*(0:(M - 1)) + pi/M;
ky = kx;
delKx = kx(2)-kx(1);
delKy = delKx;
energy = linspace(firstEnergy, lastEnergy, nEnergyPoints);
kSpaceEigenValues = zeros(M, M, 2*nOrbitals);
kSpaceEigenVectors = zeros(M, M, 2*nOrbitals, 2*nOrbitals);
kSpaceEigenValuesNormal = zeros(M, M, nOrbitals);
kSpaceEigenVectorsNormal = zeros(M, M, nOrbitals, nOrbitals);
for iKx = 1:M
        for iKy = 1:M
            k = [kx(iKx) ky(iKy)];
            % diagonalizing for normal state DOS
            kSpaceHopping = 0;
            for iUnitCell = 1:nUnitCells
                iLatticeVector = latticeVectors(iUnitCell,:);
                kSpaceHopping = kSpaceHopping + TBparameters(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
            end                        
            [eigVectorNormal eigValueNormal] = eig(kSpaceHopping);
            [eigValueKNormal sortingIndexNormal] = sort(real(diag(eigValueNormal)));
            eigVectorKNormal = (eigVectorNormal(:,sortingIndexNormal))';
            kSpaceEigenValuesNormal(iKx, iKy, :) = eigValueKNormal;
            kSpaceEigenVectorsNormal(iKx, iKy, :,:) =  eigVectorKNormal;
            % diagonalizing for SC state DOS
            kSpaceGap = 0;
            for iUnitCellDelta = 1:nUnitCellsDelta
                iLatticeVectorDelta = latticeVectorsSC(iUnitCellDelta,:);
                kSpaceGap = kSpaceGap + delta(:,:,iUnitCellDelta)*exp(1i*(iLatticeVectorDelta*k'));
            end
            kSpaceHamiltonian = [kSpaceHopping -kSpaceGap; -kSpaceGap' -kSpaceHopping];
            [eigVector eigValue] = eig(kSpaceHamiltonian);
            [eigValueK sortingIndex] = sort(real(diag(eigValue)));
            eigVectorK = (eigVector(:,sortingIndex))';
            kSpaceEigenValues(iKx, iKy, :) = eigValueK;
            kSpaceEigenVectors(iKx, iKy, :,:) =  eigVectorK;
        end
        disp(iKx)
end

% Normal state DOS
disp('Computing normal state DOS......')
greensDiagonalNormal = zeros(nOrbitals,nEnergyPoints);
countLoop = 0;
for iEnergyPoint = 1:nEnergyPoints
    for jBand = 1:nOrbitals
        greensKSpaceNormal = 0;
        E = energy(iEnergyPoint);
        for iBand = 1:nOrbitals
            % find eigenvector elements in this band
            En = squeeze(kSpaceEigenValuesNormal(:,:,iBand));
            bn = squeeze(kSpaceEigenVectorsNormal(:,:,iBand,jBand));
            % find contribution of this band to green's function 
            greensKSpaceIBand = ((abs(bn)).^2)./(E + 1i*ita - En);
            greensKSpaceNormal = greensKSpaceNormal + greensKSpaceIBand;
        end
        greensDiagonalNormal(jBand, iEnergyPoint) = (1/(2*pi))^2*delKx*delKy*singular_double_quad(1./greensKSpaceNormal);
    end
    countLoop = countLoop + 1;
    disp(countLoop);
end
bandDOSNormal = -(1/pi)*imag(greensDiagonalNormal);
totalDOSNormal = sum(bandDOSNormal);


% SC state dos
disp('Computing SC state DOS......')
eigValuesPlus = kSpaceEigenValues(:,:,(nOrbitals + 1):end); % choose positive branch of spectrum
eigVectorsPlus = kSpaceEigenVectors(:,:,(nOrbitals + 1):end,:); % corresponding eigenvectors
u = eigVectorsPlus(:,:,:,1:nOrbitals);
v = eigVectorsPlus(:,:,:,(nOrbitals+1):end);
countLoop = 0;
greensDiagonal = zeros(nOrbitals,nEnergyPoints); 
for iEnergyPoint = 1:nEnergyPoints
    for jBand = 1:nOrbitals
        greensKSpace = 0;
        E = energy(iEnergyPoint);
        for iBand = 1:nOrbitals
            % find eigenvector elements in this band
            En = squeeze(eigValuesPlus(:,:,iBand));
            un = squeeze(u(:,:,iBand,jBand));
            vn = squeeze(v(:,:,iBand,jBand));
            % find contribution of this band to green's function 
            greensKSpaceIBand = ((abs(un)).^2)./(E + 1i*ita - En) + ((abs(vn)).^2)./(E + 1i*ita + En);
            greensKSpace = greensKSpace + greensKSpaceIBand;
        end
        greensDiagonal(jBand, iEnergyPoint) = (1/(2*pi))^2*delKx*delKy*singular_double_quad(1./greensKSpace);
    end
        countLoop = countLoop + 1;
        disp(countLoop);
end
bandDOS = -(1/pi)*imag(greensDiagonal);
totalDOS = sum(bandDOS);

save DOS_homogeneous_FeSe_Milan_N_9_GammaCut_2 energy bandDOSNormal bandDOS 
% Plotting
figure; 
plot(energy,(5/nOrbitals)*totalDOSNormal,'k'); hold; plot(energy,(5/nOrbitals)*totalDOS, 'r');
axis('square'); title('Normal Vs SC dos')
figure;
plot(energy, (5/nOrbitals)*totalDOS, 'k');
hold
plot(energy, bandDOS(1,:), 'r');
plot(energy, bandDOS(2,:), 'g');
plot(energy, bandDOS(3,:), 'c');
plot(energy, bandDOS(4,:), 'm');
plot(energy, bandDOS(5,:), 'b');
axis('square'); title('Orbital resolved SC dos')