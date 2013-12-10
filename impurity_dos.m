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
    load(TB_file);
    % possibly not necessary?
    latticeVectors = latticeVector;
    % save variables for next run
    save('standart_input_imp_dos.mat')
else
    % set all the values from the inputfile
    load(inputfile);    
end;

load(TB_file);
load(Gamma_file);
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
% variable not used ?
nUnitCellsDelta = size(superDeltaVectors,1);
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
    division=str2num(division);
    part=str2num(part);
    if part > division
        % just do the summation and integration
        startindex=1;
        endindex=0;
    else
        pointspertask=ceil(M^2/division);
        startindex=pointspertask*(part-1)+1;
        endindex=startindex+pointspertask-1;
        if endindex>M^2
            % if division is not divisor of M^2, the last task has to do less
            endindex=M^2
        end;
    end
end;

% only allocate this variable if it is really needed
if (division==0 || part>division)
    greensKSpace = zeros(M, M, nDosSites, nEnergyPoints);
end;

LDOSfileName = [casestring,'_Vimp_', num2str(Vimp),  '_N_', num2str(N), '_M_', num2str(M), '_ita_', num2str(ita)];
% create some sub-directory to avoid many files in one directory
if division>0
    dirstring=[LDOSfileName,'_division_',num2str(division)];
    if ~(exist(dirstring,'dir'))
        mkdir(dirstring)
    end;
end;

% only one for loop
for index=startindex:endindex
    iKy= mod(index-1,M)+1;
    iKx= ceil(index/M);
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
    Ek_vector=eigValueK((nBands + 1):end);
    uK = eigVectorK(siteIndices,(nBands + 1):end);
    vK = eigVectorK(nBands + siteIndices,(nBands + 1):end);
    % depending on the mode do different things
    if division==0
        % single calculation of full DOS
        Ek = repmat(Ek_vector, 1, nEnergyPoints) ;
        greensKSpace(iKx, iKy, :, :) = ((abs(uK)).^2)*(1./(E - Ek + 1i*ita )) + ...
            ((abs(vK)).^2)*(1./(E + Ek + 1i*ita ));
    else
        % save result in one single file
        save([dirstring,'/','index_',num2str(index),'.mat'],'uK','vK','Ek_vector');
    end;
    toc;
end
LDOSfileName = [casestring,'_Vimp_', num2str(Vimp),  '_N_', num2str(N), '_M_', num2str(M), '_ita_', num2str(ita)];
if part>division
    if division>0
        % read in the precalculated results and sum over
        %for partindex=1:division
        %    Greenskspacefilename=[LDOSfileName,'_division_',num2str(division),'_part_',num2str(partindex)];
        %    greensKSpace_partial=load(Greenskspacefilename,'-mat');
        %    greensKSpace=greensKSpace+greensKSpace_partial.greensKSpace;
        %end;
        for index=1:M^2
            iKy= mod(index-1,M)+1;
            iKx= ceil(index/M);
            load([dirstring,'/','index_',num2str(index),'.mat']);
            % caeful: double code here, change both when doing any
            % modifications
            Ek = repmat(Ek_vector, 1, nEnergyPoints) ;
            greensKSpace(iKx, iKy, :, :) = ((abs(uK)).^2)*(1./(E - Ek + 1i*ita )) + ...
                ((abs(vK)).^2)*(1./(E + Ek + 1i*ita ));
        end;
    end
    % do the calculation of dos
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
    % to be done: change filename to general string
    save(LDOSfileName, 'energy', 'orbitalLDOSFarAway', 'orbitalLDOSImp', 'orbitalLDOSImpNN', 'orbitalLDOSImpNNN');
else
    % write out partial result for later use
    % to be done: write out only v_k, u_k
  %  Greenskspacefilename=[LDOSfileName,'_division_',num2str(division),'_part_',num2str(part)];
  %  save(Greenskspacefilename,'greensKSpace');
end
r=1;
