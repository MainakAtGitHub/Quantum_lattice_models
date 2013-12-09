% load relevant files
N = input('Enter N   ');
Vimp = input('enter impurity potential    ');
M = input('Enter no of k-points   ');% no of K points in x
ita = input('Enter ita   '); % broadening
firstEnergy = input('Enter starting energy   ');
lastEnergy = input('Enter last energy   ');
nEnergyPoints = input('Enter no of enery points   ');

load TB_hamiltonian_FeSe_2D.mat
latticeVectors = latticeVector;
load Gamma_FeSe_10_orbital_Milan_symmetrized.mat
BdGfileName = ['BdG_Impurity_FeSe', '_N_', num2str(N),'_Vimp_', num2str(Vimp)];
load(BdGfileName);


nOrbitals = size(TBparameters,1);
nBands = N^2*nOrbitals;
kx = (2*pi/M)*(0:(M - 1)) + pi/M;
ky = kx;
delKx = kx(2)-kx(1);
delKy = delKx;
energy = linspace(firstEnergy, lastEnergy, nEnergyPoints);

% lattice range for impurity, nn, nnn and far away sites
farAwayCell = [1 1];
impCell = [ceil(N/2) ceil(N/2)];
impNNCell = impCell + [0 1];
farAwaySiteIndex = ((farAwayCell(1)-1)*N + farAwayCell(2) - 1)*nOrbitals + (1:nOrbitals/2);
impSiteIndex = ((impCell(1)-1)*N + impCell(2) - 1)*nOrbitals + (1:nOrbitals/2);
impNNSiteIndex = ((impCell(1)-1)*N + impCell(2) - 1)*nOrbitals + ((nOrbitals/2+1):nOrbitals);
impNNNSiteIndex = ((impNNCell(1)-1)*N + impNNCell(2) - 1)*nOrbitals + (1:nOrbitals/2);
siteIndices = [farAwaySiteIndex impSiteIndex impNNSiteIndex impNNNSiteIndex];
nDosSites = length(siteIndices);

% Supercell quantities
maxHop = max(max(abs(latticeVectorsSC)));
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) = ...
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) - mu*eye(nOrbitals); 
[HSuper, superLatticeVectors] = supercell_hoppings(N, TBparameters, latticeVectors);
[deltaSuper,superDeltaVectors] = supercell_delta(nOrbitals, delta, maxHop);
HImpurity = zeros(nBands);
[iRange, jRange] = find_lattice_translation_index(N, nOrbitals, impCell, impCell);
impPotential = Vimp*eye(nOrbitals/2, nOrbitals/2);
HImpurity(iRange, jRange) = [impPotential zeros(nOrbitals/2); zeros(nOrbitals/2) zeros(nOrbitals/2)];

% supercell diagonalization
nSuperCells = size(superLatticeVectors,1);
nUnitCellsDelta = size(superDeltaVectors,1);
E = repmat(energy,nBands,1);
greensKSpace = zeros(M, M, nDosSites, nEnergyPoints);
for iKx = 1:M
        for iKy = 1:M
            tic;
            disp([iKx iKy]);
            k = [kx(iKx) ky(iKy)];
            kSpaceHopping = 0;
            kSpaceGap = 0;
            for iUnitCell = 1:nSuperCells
                iLatticeVector = superLatticeVectors(iUnitCell,:);
                kSpaceHopping = kSpaceHopping + HSuper(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
                kSpaceGap = kSpaceGap + deltaSuper(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
            end         
            KESuper = kSpaceHopping + HImpurity;
            kSpaceHamiltonian = [KESuper -kSpaceGap; -kSpaceGap' -KESuper];
            [eigVector eigValue] = eig(kSpaceHamiltonian);
            [eigValueK sortingIndex] = sort(real(diag(eigValue)));
            eigVectorK = (eigVector(:,sortingIndex));
            Ek = repmat(eigValueK((nBands + 1):end), 1, nEnergyPoints) ;
            uK = eigVectorK(siteIndices,(nBands + 1):end);
            vK = eigVectorK(nBands + siteIndices,(nBands + 1):end);
            greensKSpace(iKx, iKy, :, :) = ((abs(uK)).^2)*(1./(E - Ek + 1i*ita )) + ...
                                           ((abs(vK)).^2)*(1./(E + Ek + 1i*ita ));
            toc;           
        end
end

greensRealSpace = zeros(nDosSites, nEnergyPoints);
for iSite = 1: nDosSites
    for iEnergyPoint = 1:nEnergyPoints
        greensRealSpace(iSite, iEnergyPoint) = (1/(2*pi))^2*delKx*delKy*...
        singular_double_quad(1./squeeze(greensKSpace(:, :, iSite, iEnergyPoint)));
    end
end

ldos = (-(1/pi))*imag(greensRealSpace);
orbitalLDOSFarAway = ldos(1:5,:);
totalLDOSFarAway = sum(orbitalLDOSFarAway,1);
orbitalLDOSImp = ldos(6:10,:);
totalLDOSImp = sum(orbitalLDOSImp,1);
orbitalLDOSImpNN = ldos(11:15,:);
totalLDOSImpNN = sum(orbitalLDOSImpNN,1);
orbitalLDOSImpNNN = ldos(16:20,:);
totalLDOSImpNNN = sum(orbitalLDOSImpNNN,1);   

LDOSfileName = ['LDOS_FeSe_Milan_Gamma','_Vimp_', num2str(Vimp),  '_N_', num2str(N), '_M_', num2str(M), '_ita_', num2str(ita)];
save(LDOSfileName, 'energy', 'orbitalLDOSFarAway', 'orbitalLDOSImp', 'orbitalLDOSImpNN', 'orbitalLDOSImpNNN');
