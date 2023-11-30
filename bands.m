function h=bands(inputfile,calcSC)

% Modified homogeneous_dos.m
% takes \Delta_ij as input and constructs \Delta_i0.
if ~exist('calcSC','var')
    calcSC=false;
end;
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
if calcSC
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
else
    mu=0;
end;


nUnitCells = size(latticeVector,1);
%TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) = ...
%TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) - mu*eye(nOrbitals);

% setup of a path
numpoints=40;
kx1=[0:pi/40:pi];
ky1=0*kx1;
ky2=[pi-pi/40:-pi/40:0];
%[kx,ky]=meshgrid(-pi:pi/20:pi,-pi:pi/20:pi);
%kx=kx(:)';
%ky=ky(:)';
kx=[kx1,pi*ones(1,40),ky2];
ky=[kx1,ky2,0*ky2];
numk=size(kx,2);
kSpaceEigenValuesNormal = zeros(numk, nOrbitals);
kSpaceEigenVectorsNormal = zeros(numk, nOrbitals, nOrbitals);
for iK = 1:numk
            k = [kx(iK) ky(iK)];
            % diagonalizing for normal state DOS
            kSpaceHopping = 0;
            for iUnitCell = 1:nUnitCells
                %convert to 2D lattice vector (assuming third dimension is
                %zero)
                iLatticeVector = latticeVector(iUnitCell,1:2);
                exp(1i*(iLatticeVector*k'));
                kSpaceHopping = kSpaceHopping + TBparameters(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
            end                        
            [eigVectorNormal eigValueNormal] = eig(kSpaceHopping);
            [eigValueKNormal sortingIndexNormal] = sort(real(diag(eigValueNormal)));
            eigVectorKNormal = (eigVectorNormal(:,sortingIndexNormal))';
            kSpaceEigenValuesNormal(iK, :) = eigValueKNormal;
            kSpaceEigenVectorsNormal(iK, :,:) =  eigVectorKNormal;
            % diagonalizing for SC state DOS
            if calcSC
            kSpaceGap = 0;
            for iUnitCellDelta = 1:nUnitCellsDelta
                iLatticeVectorDelta = latticeVectorsSC(iUnitCellDelta,:);
                kSpaceGap = kSpaceGap + delta(:,:,iUnitCellDelta)*exp(1i*(iLatticeVectorDelta*k'));
            end
            kSpaceHamiltonian = [kSpaceHopping -kSpaceGap; -kSpaceGap' -kSpaceHopping];
            [eigVector eigValue] = eig(kSpaceHamiltonian);
            [eigValueK sortingIndex] = sort(real(diag(eigValue)));
            eigVectorK = (eigVector(:,sortingIndex))';
            kSpaceEigenValues(iK, :) = eigValueK;
            kSpaceEigenVectors(iK, :,:) =  eigVectorK;
            end;
        if  mod(iK,10)==0   
            disp(['Done ',num2str(iK), ' of ',num2str(M),' kx values.']);
        end;
end
figure;
plot(kSpaceEigenValuesNormal);
return
save(LDOSfileName, 'energy', 'bandDOSNormal', 'bandDOS', '-mat');
if (usejava('jvm') && ~feature('ShowFigureWindows'))
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
end;