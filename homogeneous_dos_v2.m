function h=homogeneous_dos_v2(inputfile)

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

if (~exist('Vimp','var'))
    if abs(Vimp)>0
        disp('Warning: finite impurity potential, not homogeneous case.')
    end;
end;

load(TB_file,'-mat');
load(Gamma_file,'-mat');
load(BdGfileName,'-mat');


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


nUnitCells = size(latticeVector,1);
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) = ...
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) - mu*eye(nOrbitals);

M=N*M;
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
                iLatticeVector = latticeVector(iUnitCell,:);
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
if (~exist('tetra','var'))
    tetra=false;
end;

% Normal state DOS
disp('Computing normal state DOS......')
greensDiagonalNormal = zeros(nOrbitals,nEnergyPoints);    
if ~tetra
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
else
            disp('...using 2D version of Tetrahedron method');
            % set up a k-mesh that is suitable to cover the whole
            % Brillouinzone with triangles
            mesh1=[0:1:(M)]*2*pi/M;
            % to do: kx,ky can be only a vector to simplify indexing
            [kx,ky] = meshgrid(mesh1, mesh1);
                % to do: vectorize the code!
                for iband=1:nOrbitals
                    disp(['Band ',num2str(iband),' of ',num2str(nOrbitals)]);
                    E=kSpaceEigenValuesNormal(:,:,iband);
                    a=squeeze(kSpaceEigenVectorsNormal(:,:,iband,:)).*conj(squeeze(kSpaceEigenVectorsNormal(:,:,iband,:)));
                    greensDiagonalNormal(:, :) = greensDiagonalNormal(:, :) + f(E,a,kx,ky,energy);
                end;
                
bandDOSNormal = greensDiagonalNormal;

end;
totalDOSNormal = sum(bandDOSNormal);


% SC state dos
disp('Computing SC state DOS......')
eigValuesPlus = kSpaceEigenValues(:,:,(nOrbitals + 1):end); % choose positive branch of spectrum
eigVectorsPlus = kSpaceEigenVectors(:,:,(nOrbitals + 1):end,:); % corresponding eigenvectors
u = eigVectorsPlus(:,:,:,1:nOrbitals);
v = eigVectorsPlus(:,:,:,(nOrbitals+1):end);
countLoop = 0;
greensDiagonal = zeros(nOrbitals,nEnergyPoints); 
if ~tetra
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

else
            disp('...using 2D version of Tetrahedron method');
            % set up a k-mesh that is suitable to cover the whole
            % Brillouinzone with triangles
            mesh1=[0:1:(M)]*2*pi/M;
            % to do: kx,ky can be only a vector to simplify indexing
            [kx,ky] = meshgrid(mesh1, mesh1);
                % to do: vectorize the code!
                for iband=1:nOrbitals
                    disp(['Band ',num2str(iband),' of ',num2str(nOrbitals)]);
		    E=squeeze(eigValuesPlus(:,:,iband));
		    un = squeeze(u(:,:,iband,:));
		    vn = squeeze(v(:,:,iband,:));
		    a=un.*conj(un);
                    greensDiagonal(:, :) = greensDiagonal(:, :) + f(E,a,kx,ky,energy);
                    a=vn.*conj(vn);
                    greensDiagonal(:, :) = greensDiagonal(:, :) + f(E,a,kx,ky,-energy);
                end;
                bandDOS=greensDiagonal;
end;
totalDOS = sum(bandDOS);

LDOSfileName0 = [casestring,'_Vimp_', num2str(Vimp),  '_N_', num2str(N)];

if ~tetra
    LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_ita_', num2str(ita)];
else
    LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_tetra_corr']
end;

disp('Writing out k-space calculated DOS ...');
save(LDOSfileName, 'energy', 'bandDOSNormal', 'bandDOS', '-mat');
if usejava('jvm') && ~feature('ShowFigureWindows')
    disp(['please plot the result using plot_homogeneous_dos_v2(''',LDOSfileName,''')']);
else
    %# GUI available
% Plotting
figure; 
plot(energy,(5/nOrbitals)*totalDOSNormal,'k'); hold; plot(energy,(5/nOrbitals)*totalDOS, 'r');
axis('square'); title('Normal Vs SC dos')
% Create legend
legend show
figure;
plot(energy, (5/nOrbitals)*totalDOS, 'k');
hold
plot(energy, bandDOS(1,:), 'r');
plot(energy, bandDOS(2,:), 'g');
plot(energy, bandDOS(3,:), 'c');
plot(energy, bandDOS(4,:), 'm');
plot(energy, bandDOS(5,:), 'b');
axis('square'); title('Orbital resolved SC dos')
% Create legend
legend show

