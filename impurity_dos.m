function r=impurity_dos(inputfile, division, part)

% calculate the impurity density of states
% inputfile : file that contains the parameters
% division : divide task into division parts
% part : calculate this part part = 1...division

% if no inputfile is specified, read in from the command line
if nargin <1
    % load relevant files
    N = input('Enter N   ');
    Vimp = input('enter impurity potential    ');
    M = input('Enter no of k-points   ');% no of K points in x
    ita = input('Enter ita   '); % broadening
    firstEnergy = input('Enter starting energy   ');
    lastEnergy = input('Enter last energy   ');
    nEnergyPoints = input('Enter no of energy points   ');
    % set some input filenames
    TB_file='TB_hamiltonian_FeSe_2D.mat';
    Gamma_file='Gamma_FeSe_10_orbital_Milan_symmetrized.mat';
    BdGfileName = ['BdG_Impurity_FeSe', '_N_', num2str(N),'_Vimp_', num2str(Vimp)];
    casestring='LDOS_FeSe_Milan_Gamma';
    % save variables for next run
    save('standart_input_imp_dos.mat')
else
    % otherwise read inputfile
    try
        % old input format with mat-file
        load(inputfile);
    catch 
        % new text-based input format
        read_input_file=inputfile;
        read_input;
        read_input_file
    end;
end;

load(TB_file);
% possibly not necessary?
%latticeVectors = latticeVector;
load(Gamma_file);
load(BdGfileName);

nOrbitals = size(TBparameters,1);
nBands = N^2*nOrbitals;
kx = (2*pi/M)*(0:(M - 1)) + pi/M;
ky = kx;
delKx = kx(2)-kx(1);
delKy = delKx;
% introduce some indexing for the k-vectors to use precalculated values
kx_ind=[0:(M-1);ones(1,M)*M];
kx_ind=kx_ind./repmat(gcd(kx_ind(1,:),kx_ind(2,:)),2,1);
ky_ind=kx_ind;
% find the
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
[HSuper, superLatticeVectors] = supercell_hoppings(N, TBparameters, latticeVector);
[deltaSuper,superDeltaVectors] = supercell_delta(nOrbitals, delta, maxHop);
HImpurity = zeros(nBands);
[iRange, jRange] = find_lattice_translation_index(N, nOrbitals, impCell, impCell);
impPotential = Vimp*eye(nOrbitals/2, nOrbitals/2);
HImpurity(iRange, jRange) = [impPotential zeros(nOrbitals/2); zeros(nOrbitals/2) zeros(nOrbitals/2)];

% supercell diagonalization
nSuperCells = size(superLatticeVectors,1);
% variable not used ?
%nUnitCellsDelta = size(superDeltaVectors,1);
E = repmat(energy,nBands,1);

% orphan later
if ~exist('casestring','var')
    casestring='LDOS_FeSe_Milan_Gamma';
end;

% some code for parallelization
if nargin < 2
    startindex=1;
    endindex=M^2;
    % do also the integration
    part=1;
    division=0;
else
    % set up the start and endindex for parallelization
    if isa(division,'char')
        division=str2num(division);
    end;
    if isa(part,'char')
        part=str2num(part);
    end
    if part > division
        % just do the summation and integration
        startindex=1;
        endindex=0;
    elseif part <0
        % only calculate single point
        startindex=-part;
        endindex=-part;
    else
        pointspertask=ceil(M^2/division);
        startindex=pointspertask*(part-1)+1;
        endindex=startindex+pointspertask-1;
        if endindex>M^2
            % if division is not divisor of M^2, the last task has to do less
            endindex=M^2
        end;
    end
end

if (~exist('tetra','var'))
    tetra=false;
end;

% only allocate this variable if it is really needed
if (division==0 || part>division)
    if ~tetra
        greensKSpace = zeros(M, M, nDosSites, nEnergyPoints);
    else
        % some huge arrays to store the result
        % store the edges twice to calculate the whole area
        ukall=zeros(M+1,M+1,nDosSites,nBands);
        vkall=zeros(M+1,M+1,nDosSites,nBands);
        ekall=zeros(M+1,M+1,nDosSites,nBands);
    end;
end;

LDOSfileName0 = [casestring,'_Vimp_', num2str(Vimp),  '_N_', num2str(N)];
% create some sub-directory to avoid many files in one directory
if division>0
    dirstring=['data_',LDOSfileName0];
    if ~(exist(dirstring,'dir'))
        mkdir(dirstring)
    end;
end;

% only one for loop
for index=startindex:endindex
    iKy= mod(index-1,M)+1;
    iKx= ceil(index/M);
    tic;
    if division > 0
        ekukvk_file=[dirstring,'/','kx_',num2str(kx_ind(1,iKx)),'_',num2str(kx_ind(2,iKx)),'ky_',num2str(ky_ind(1,iKy)),'_',num2str(ky_ind(2,iKy)),'.mat'];
    end;
    if (~exist(ekukvk_file, 'file') || division==0)
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
        clear eigValue
        eigVectorK = (eigVector(:,sortingIndex));
        Ek_vector=eigValueK((nBands + 1):end);
        clear eigValueK
        uK = eigVectorK(siteIndices,(nBands + 1):end);
        vK = eigVectorK(nBands + siteIndices,(nBands + 1):end);
        clear eigVectorK
        % depending on the mode do different things
        if division==0
            % single calculation of full DOS
            Ek = repmat(Ek_vector, 1, nEnergyPoints) ;
            greensKSpace(iKx, iKy, :, :) = ((abs(uK)).^2)*(1./(E - Ek + 1i*ita )) + ...
                ((abs(vK)).^2)*(1./(E + Ek + 1i*ita ));
        else
            % save result in one single file
            % put k_vector in filename to avoid double calculation ?
            save(ekukvk_file,'uK','vK','Ek_vector');
        end;
    else
        disp([ekukvk_file,' already calculated, skipping.' ])
    end
    toc;
end
if ~tetra
    LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_ita_', num2str(ita)];
else
    LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_tetra']
end;
if part>division
    if division>0
        % read in the precalculated results and sum over
        %for partindex=1:division
        %    Greenskspacefilename=[LDOSfileName,'_division_',num2str(division),'_part_',num2str(partindex)];
        %    greensKSpace_partial=load(Greenskspacefilename,'-mat');
        %    greensKSpace=greensKSpace+greensKSpace_partial.greensKSpace;
        %end;
        disp('Reading in precalculated eigenvalues and Bogoliubov coefficients...')
        for index=1:M^2
            iKy= mod(index-1,M)+1;
            iKx= ceil(index/M);
            if iKy==1
                disp(['reading k-point',num2str(iKx),' ',num2str(iKy)])
            end;
            try
                ekukvk_file=[dirstring,'/','kx_',num2str(kx_ind(1,iKx)),'_',num2str(kx_ind(2,iKx)),'ky_',num2str(ky_ind(1,iKy)),'_',num2str(ky_ind(2,iKy)),'.mat'];
                load(ekukvk_file); 
            catch exception
                % missing k-point (or wrong input as number of k-points)
                % First can happen if one job crashes; catch this by
                % calculating on the fly
                disp('Missing k-point, recalculating on the fly.')
                impurity_dos(inputfile, division, -index);
                load(ekukvk_file);
            end
            % caeful: double code here, change both when doing any
            % modifications
            if ~tetra
                Ek = repmat(Ek_vector, 1, nEnergyPoints) ;
                greensKSpace(iKx, iKy, :, :) = ((abs(uK)).^2)*(1./(E - Ek + 1i*ita )) + ...
                    ((abs(vK)).^2)*(1./(E + Ek + 1i*ita ));
            else
                ukall(iKx,iKy,:,:)=uK;
                vkall(iKx,iKy,:,:)=vK;
                Ekall(iKx,iKy,:)=Ek_vector;
                % store the edges twice to construct set of triangles that cover the whole area
                if iKx==1
                   ukall(M+1,iKy,:,:)=uK;
                   vkall(M+1,iKy,:,:)=vK;
                   eKALL(M+1,iKy,:,:)=Ek_vector;
                end;
                if iKy==1
                   ukall(iKx,M+1,:,:)=uK;
                   vkall(iKx,M+1,:,:)=vK;
                   eKALL(iKx,M+1,:,:)=Ek_vector;
                end;
                if (iKx==1) && (iKy==1)
                   ukall(M+1,M+1,:,:)=uK;
                   vkall(M+1,M+1,:,:)=vK;
                   eKALL(M+1,M+1,:,:)=Ek_vector;
                end;
            end;
        end;
        disp('... done.');
    end
    % do the calculation of dos
    greensRealSpace = zeros(nDosSites, nEnergyPoints);
    disp('calculating GF in real space...');
    for iSite = 1: nDosSites
        disp(['Done ',num2str(iSite),' of ', num2str(nDosSites), 'nDosSites']);
        if ~tetra
            for iEnergyPoint = 1:nEnergyPoints            
                greensRealSpace(iSite, iEnergyPoint) = (1/(2*pi))^2*delKx*delKy*...
                    singular_double_quad(1./squeeze(greensKSpace(:, :, iSite, iEnergyPoint)));
            end
        else
            disp('...using 2D version of Tetrahedron method');
            mesh1=[0.5:1:(M-0.5)]*2*pi/M;
            [kx,ky] = meshgrid(mesh1, mesh1);
                % to do: vectorize the code!
                for iband=1:nBands
                    disp(['Band ',num2str(iband),' of ',num2str(nBands)]);
                    E=Ekall(:,:,iband);
                    a=ukall(:,:,iSite,iband).*conj(ukall(:,:,iSite,iband));
                    greensRealSpace(iSite, :) = greensRealSpace(iSite, :) + f(E,a,kx,ky,energy);
                    a=vkall(:,:,iSite,iband).*conj(vkall(:,:,iSite,iband));
                    greensRealSpace(iSite, :) = greensRealSpace(iSite, :) + f(E,a,kx,ky,-energy);
                end;
        end;
    end
    disp('Writing out LDOS ...');
    if tetra
        ldos=greensRealSpace;
    else
        ldos = (-(1/pi))*imag(greensRealSpace);
    end;
    orbitalLDOSFarAway = ldos(1:5,:);
    totalLDOSFarAway = sum(orbitalLDOSFarAway,1); %#ok<NASGU>
    orbitalLDOSImp = ldos(6:10,:);
    totalLDOSImp = sum(orbitalLDOSImp,1);%#ok<NASGU>
    orbitalLDOSImpNN = ldos(11:15,:);
    totalLDOSImpNN = sum(orbitalLDOSImpNN,1);%#ok<NASGU>
    orbitalLDOSImpNNN = ldos(16:20,:);
    totalLDOSImpNNN = sum(orbitalLDOSImpNNN,1);%#ok<NASGU>
    % to be done: change filename to general string
    save(LDOSfileName, 'energy', 'orbitalLDOSFarAway', 'orbitalLDOSImp', 'orbitalLDOSImpNN', 'orbitalLDOSImpNNN');
end
r=1;
