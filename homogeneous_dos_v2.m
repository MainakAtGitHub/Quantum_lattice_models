function h=homogeneous_dos_v2(inputfile,calcSC)

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
if ~exist('calcSC','var')
    calcSC=true;
end;
if (~exist('Vimp','var'))
    if abs(Vimp)>0
        disp('Warning: finite impurity potential, not homogeneous case.')
    end;
end;
if ~(exist('singular_quad','var'))
    singular_quad=true;
end;
if ~(exist('write_states','var'))
    write_states=false;
end;
if ~(exist('calc_dos','var'))
    calc_dos=true;
end;
if ~singular_quad
        sqstring='sum';
else
    sqstring='';
end;
%%%%Mainak
if ~exist('pos_file','var')
    load(TB_file,'-mat');
    nOrbitals = size(TBparameters,1);
else
    if exist('further_N_cut_off','var')
        %%%%%%%Mainak
        if further_N_cut_off > 2
            latticeVector=[1,0,0;...
                -1,0,0;...
                0,1,0;...
                0,-1,0;...
                1,1,0;...
                1,-1,0;...
                -1,1,0;...
                -1,-1,0;...
                -2,0,0;...
                2,0,0;...
                0,-2,0;...
                0,2,0;...
                0,0,0];   %FOR NN and NNN and NNNN
        elseif and(further_N_cut_off > sqrt(2),further_N_cut_off < 2)
            latticeVector=[1,0,0;...
                -1,0,0;...
                0,1,0;...
                0,-1,0;...
                1,1,0;...
                1,-1,0;...
                -1,1,0;...
                -1,-1,0;...
                0,0,0];   %FOR NN and NNN
        else
            latticeVector=[1,0,0;-1,0,0;0,1,0;0,-1,0;0,0,0];
        end
    else
        latticeVector=[1,0,0;-1,0,0;0,1,0;0,-1,0;0,0,0];
    end
    t_r_ref=load(ref_grid_hopping_file);
    TBparameters=zeros(1,1,size(latticeVector,1));
    for itr_TB = 1:size(latticeVector,1)
        if exist('ref_grid_hopping_file','var')
            TBparameters(:,:,itr_TB)=griddata(t_r_ref(:,1),t_r_ref(:,2),10^-3*t_r_ref(:,3),latticeVector(itr_TB,1),latticeVector(itr_TB,2),'natural');
        else
            %%%% TB parameters for the toy hopping, gaussian in distance (below)
            TBparameters(:,:,itr_TB)=-1/exp(-1)*exp(-((latticeVector(itr_TB,1))^2+(latticeVector(itr_TB,2))^2));
        end
        if itr_TB == size(latticeVector,1)
            TBparameters(:,:,itr_TB)=0;
        end
    end
    nOrbitals = size(TBparameters,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Mainak
if calcSC
    load(Gamma_file,'-mat');
    load(BdGfileName,'-mat');
    N = sqrt(size(delta,1)/nOrbitals);
    if exist('latticeVectorsSC','var')
        nUnitCellsDelta = size(latticeVectorsSC,1);
    else
        latticeVectorsSC=Gammafull.latt;
        nUnitCellsDelta = size(Gammafull.latt,1);
    end
    deltaCenter = zeros(nOrbitals, nOrbitals, nUnitCellsDelta);
    jCell = [ceil(N/2) ceil(N/2)];
    for i = 1:nUnitCellsDelta
        iCell = jCell + latticeVectorsSC(i,:);
        [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, iCell, jCell);
        deltaCenter(:,:,i) = delta(iRange, jRange);
    end
    delta = deltaCenter;
else
    load(BdGfileName,'mu','-mat');
%   mu=0;
end;


% fix for missing on-site energy in tb parameters below
%if sum((latticeVector(:,1)==0) & (latticeVector(:,2)==0))==0
 %   TBparameters(:,:,end+1)=- mu*eye(nOrbitals);
 %   latticeVector(end+1,:)=0*latticeVector(1,:);
%end;
nUnitCells = size(latticeVector,1);
%TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) = ...
%TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) - mu*eye(nOrbitals);

M=N*M;
kx = (2*pi/M)*(0:(M - 1));% + pi/M;
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
            kSpaceHopping = - mu*eye(nOrbitals);
            kSpaceHoppingc = - mu*eye(nOrbitals);
            %zeros(nOrbitals,nOrbitals);
            for iUnitCell = 1:nUnitCells
                iLatticeVector = latticeVector(iUnitCell,1:2);
                kSpaceHopping = kSpaceHopping + TBparameters(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
                kSpaceHoppingc = kSpaceHoppingc + TBparameters(:,:,iUnitCell)*exp(-1i*(iLatticeVector*k'));
            end                        
            [eigVectorNormal, eigValueNormal] = eig(kSpaceHopping);
            [eigValueKNormal, sortingIndexNormal] = sort(real(diag(eigValueNormal)));
            eigVectorKNormal = (eigVectorNormal(:,sortingIndexNormal))';
            kSpaceEigenValuesNormal(iKx, iKy, :) = eigValueKNormal;
            kSpaceEigenVectorsNormal(iKx, iKy, :,:) =  eigVectorKNormal;
            % diagonalizing for SC state DOS
            if calcSC
            kSpaceGap = zeros(nOrbitals,nOrbitals);
            for iUnitCellDelta = 1:nUnitCellsDelta
                iLatticeVectorDelta = latticeVectorsSC(iUnitCellDelta,:);
                kSpaceGap = kSpaceGap + delta(:,:,iUnitCellDelta)*exp(1i*(iLatticeVectorDelta*k'));
            end
            kSpaceHamiltonian = [kSpaceHopping -kSpaceGap; -kSpaceGap' -conj(kSpaceHoppingc)];
            [eigVector, eigValue] = eig(kSpaceHamiltonian);
            [eigValueK, sortingIndex] = sort(real(diag(eigValue)));
            eigVectorK = (eigVector(:,sortingIndex))';
            kSpaceEigenValues(iKx, iKy, :) = eigValueK;
            kSpaceEigenVectors(iKx, iKy, :,:) =  eigVectorK;
            end;
        end
        if  mod(iKx,10)==0   
            disp(['Done ',num2str(iKx), ' of ',num2str(M),' kx values.']);
        end;
end
if (~exist('tetra','var'))
    tetra=false;
end;
% save the eigenvalues for later processing
if write_states
    outfile1=[BdGfileName,'_M_', num2str(M)];
    save([outfile1,'_normal'], 'kSpaceEigenValuesNormal', 'kx', 'ky', '-mat');
    save([outfile1,'_SC'], 'kSpaceEigenValues', 'kx', 'ky', '-mat');
end;
if calc_dos
% Normal state DOS
disp('Computing normal state DOS......')
greensDiagonalNormal = zeros(nOrbitals,nEnergyPoints);    
if ~tetra
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
        if singular_quad 
            %greensDiagonalNormal(jBand, iEnergyPoint) = (1/(2*pi))^2*delKx*delKy*singular_double_quad(1./greensKSpaceNormal);
            greensDiagonalNormal(jBand, iEnergyPoint) = (1/(2*pi))^2*delKx*delKy*singular_double_quad_mex(add_BZ_boundary(1./greensKSpaceNormal));
        else
            greensDiagonalNormal(jBand, iEnergyPoint) = (1/(2*pi))^2*delKx*delKy*sum(sum(greensKSpaceNormal));
        end;
    end
    if  mod(iEnergyPoint,10)==0
        disp(['Done ',num2str(iEnergyPoint), ' of ',num2str(nEnergyPoints)]);
    end;
end
bandDOSNormal = -(1/pi)*imag(greensDiagonalNormal);
else
            disp('...using 2D version of Tetrahedron method');
            % set up a k-mesh that is suitable to cover the whole
            % Brillouinzone with triangles
            mesh1=[0:1:(M)+.5]*2*pi/M;
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
totalDOSNormal = sum(bandDOSNormal,1);

if calcSC
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
        if singular_quad
           %greensDiagonal(jBand, iEnergyPoint) = (1/(2*pi))^2*delKx*delKy*singular_double_quad(1./greensKSpace);
            greensDiagonal(jBand, iEnergyPoint) = (1/(2*pi))^2*delKx*delKy*singular_double_quad_mex(add_BZ_boundary(1./greensKSpace));
        else
            greensDiagonal(jBand, iEnergyPoint) = (1/(2*pi))^2*delKx*delKy*sum(sum(greensKSpace));
        end;
    end
        countLoop = countLoop + 1; 
        if  mod(iEnergyPoint,10)==0
            disp(['Done ',num2str(countLoop), ' of ',num2str(nEnergyPoints)]);
        end;
end
bandDOS = -(1/pi)*imag(greensDiagonal);

else
            disp('...using 2D version of Tetrahedron method');
            % set up a k-mesh that is suitable to cover the whole
            % Brillouinzone with triangles
            mesh1=[0:1:(M)+.5]*2*pi/M;
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
totalDOS = sum(bandDOS,1);

LDOSfileName0 = [casestring,'_Vimp_', num2str(Vimp),  '_N_', num2str(N)];
else
    LDOSfileName0 = [casestring,'normal_Vimp_', num2str(Vimp),  '_N_', num2str(N)];
    bandDOS=[];
end;
if ~tetra
    LDOSfileName = [LDOSfileName0,sqstring, '_M_', num2str(M),'_ita_', num2str(ita)];
else
    LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_tetra_corr1']
end;

disp('Writing out k-space calculated DOS ...');
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
%%Mainak
% legend show
%%Mainak
figure;
plot(energy, (5/nOrbitals)*totalDOS, 'k');
hold
plot(energy, bandDOS(1,:), 'r');
try
plot(energy, bandDOS(2,:), 'g');
plot(energy, bandDOS(3,:), 'c');
plot(energy, bandDOS(4,:), 'm');
plot(energy, bandDOS(5,:), 'b');
catch
end;
axis('square'); title('Orbital resolved SC dos')
% Create legend
%%Mainak
% legend show
%%Mainak
end;
end;

