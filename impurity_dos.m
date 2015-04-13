function r=impurity_dos(inputfile, division, part)

% calculate the impurity density of states
% inputfile : file that contains the parameters
% division : divide task into division parts
% part : calculate this part part = 1...division
calcGreens =true; % set default to calculation of Greens function
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
sqstring='';
if ~(exist('singular_quad','var'))
    singular_quad=true;
end;
if ~(exist('sublattice','var'))
    sublattice=1;
end;
if ~singular_quad
        sqstring='sum';
end;
if ~(exist('debug','var'))
    debug=false;
end;
if ~(exist('nDosSitesfile','var'))
    nDosSitesfile='';
end;
if ~(exist('version','var'))
    version=1.0;
end;
[dirprefix,~,~] = fileparts(inputfile);
if ~isempty(dirprefix)
    dirprefix=[dirprefix,filesep]
end;
load(TB_file);
% possibly not necessary?
%latticeVectors = latticeVector;
load(Gamma_file);
% explicitely overwrite function mu.m from Matlab
mu=0;
load(BdGfileName);
nOrbitals = size(TBparameters,1);
nBands = N^2*nOrbitals;
% for reuse of the k-points with different M, we dropp the shift of the k-grid!
kx = (2*pi/M)*(0:(M - 1));
ky = kx;
delKx = 2*pi/M; %kx(2)-kx(1);
delKy = delKx;
% introduce some indexing for the k-vectors to use precalculated values
kx_ind=[0:(M-1);ones(1,M)*M];
kx_ind=kx_ind./repmat(gcd(kx_ind(1,:),kx_ind(2,:)),2,1);
ky_ind=kx_ind;


if ~calcGreens
    if isempty(nDosSitesfile)
        % lattice range for impurity, nn, nnn and far away sites
        farAwayCell = [1 1];
        impCell = [ceil(N/2) ceil(N/2)];
        impNNCell = impCell + [0 1];
        switch sublattice
            case 1
                efforb=nOrbitals/2;
                farAwaySiteIndex = ((farAwayCell(1)-1)*N + farAwayCell(2) - 1)*nOrbitals + (1:nOrbitals/2);
                impSiteIndex = ((impCell(1)-1)*N + impCell(2) - 1)*nOrbitals + (1:nOrbitals/2);
                impNNSiteIndex = ((impCell(1)-1)*N + impCell(2) - 1)*nOrbitals + ((nOrbitals/2+1):nOrbitals);
                impNNNSiteIndex = ((impNNCell(1)-1)*N + impNNCell(2) - 1)*nOrbitals + (1:nOrbitals/2);
            case -1
                % not fully correct, but works fine for C_4 symmetric
                % impurity
                efforb=nOrbitals/2;
                farAwaySiteIndex = ((farAwayCell(1)-1)*N + farAwayCell(2) - 1)*nOrbitals + (1:nOrbitals/2);
% upper orbitals are now on impurity
                impSiteIndex =  ((impCell(1)-1)*N + impCell(2) - 1)*nOrbitals + ((nOrbitals/2+1):nOrbitals);
% lower orbitals are on the NN site
                impNNSiteIndex=((impCell(1)-1)*N + impCell(2) - 1)*nOrbitals + (1:nOrbitals/2);
                %impNNSiteIndex = ((impCell(1)-1)*N + impCell(2) - 1)*nOrbitals + ((nOrbitals/2+1):nOrbitals);
                impNNNSiteIndex = ((impNNCell(1)-1)*N + impNNCell(2) - 1)*nOrbitals +  ((nOrbitals/2+1):nOrbitals);
            case 0
                efforb=nOrbitals;
                farAwaySiteIndex = ((farAwayCell(1)-1)*N + farAwayCell(2) - 1)*nOrbitals + (1:nOrbitals);
                impSiteIndex = ((impCell(1)-1)*N + impCell(2) - 1)*nOrbitals + (1:nOrbitals);
                impNNSiteIndex = ((impNNCell(1)-1)*N + impNNCell(2) - 1)*nOrbitals + (1:nOrbitals);
                impNNNCell = impCell + [1 1];
                impNNNSiteIndex = ((impNNNCell(1)-1)*N + impNNNCell(2) - 1)*nOrbitals + (1:nOrbitals);
        end;
        siteIndices = [farAwaySiteIndex impSiteIndex impNNSiteIndex impNNNSiteIndex];
        LDOSsites=[-inf,-inf; 0,0; 0,1;1,1];
        nDosSites = length(siteIndices);
    else
        % user defined sites
        load(nDosSitesfile);
        % set up the siteIndices
        siteIndices=[];
        impCell = [ceil(N/2) ceil(N/2)];
        for n=1:size(LDOSsites,1)
            if abs(sublattice)>0
                efforb=nOrbitals/2;
                disp('not fully tested for 10 orbitals yet, works modulo C_4 symmetry')
                a=sum(LDOSsites(n,:))/2;
                b=LDOSsites(n,2)-LDOSsites(n,2);
                cell2=[ceil(b) floor(a)];
                if LDOSsites(n,1)==-inf
                    cell=[ 1 1];
                else
                    cell=impCell+cell2;
                end
                if floor(a)==a
                    siteIndex=((cell(1)-1)*N + cell(2) - 1)*nOrbitals +(1:efforb);
                else
                    siteIndex=((cell(1)-1)*N + cell(2) - 1)*nOrbitals +(1:efforb)+efforb;
                end;
            else
                efforb=nOrbitals;
                if LDOSsites(n,1)==-inf
                    % far away point!
                    cell=[ 1 1];
                else
                    cell=impCell+LDOSsites(n,:);
                end
                siteIndex=((cell(1)-1)*N + cell(2) - 1)*nOrbitals +(1:efforb);
            end
            siteIndices=[siteIndices siteIndex];
        end
        nDosSites = length(siteIndices);
    end;
else
    % calcGreens=true needs tetra=false
    tetra=false
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

% only set up quantities if needed
if ~(division==0 || part>division)
% Supercell quantities
if exist('latticeVectorsSC','var')
	maxHop = max(max(abs(latticeVectorsSC)));
else
	maxHop = max(max(abs(Gammafull.latt)));
end
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) = ...
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) - mu*eye(nOrbitals); 
[HSuper, superLatticeVectors] = supercell_hoppings(N, TBparameters, latticeVector);
% similar code as in BdG_impurity to make the hoppings numerically a
% Hermitean matrix
HSuper=0.5*(HSuper+HSuper');
[deltaSuper,superDeltaVectors] = supercell_delta(nOrbitals, delta, maxHop);
%deltaSuper=sparse(deltaSuper);
%HSuper=sparse(HSuper);
% to be done: implementation of more complicated impurity potentials

% supercell diagonalization
nSuperCells = size(superLatticeVectors,1);
end
% variable not used ?
%nUnitCellsDelta = size(superDeltaVectors,1);
if ~calcGreens
    % find the
energy = linspace(firstEnergy, lastEnergy, nEnergyPoints);
if version >0
    E = repmat(energy,nBands*2,1);
else
E = repmat(energy,nBands,1);
end
% only one loop below
le=1;
else
    if ~(exist('Greensenergy','var'))
        disp('No Greensenergy given, setting to 0.');
        E=0;
    else
        if ~ischar(Greensenergy)
            E=Greensenergy;
            le=1;
        else
            % load the list of energies to be calculated from the given
            % file
            load(Greensenergy,'E');
            le=numel(E);
        end;
    end
end;

% orphan later
if ~exist('casestring','var')
    casestring='LDOS_FeSe_Milan_Gamma';
end;

if (~exist('tetra','var'))
    tetra=false;
end;

% only allocate this variable if it is really needed
if (division==0 || part>division)
    if ~tetra
        if ~calcGreens
            greensKSpace = zeros(M, M, nDosSites, nEnergyPoints);
        else
            latticeGreensK = zeros(M, M, nBands, nBands);
        end
    else
        % some huge arrays to store the result
        % store the edges twice to calculate the whole area
        if version >0
        ukall=zeros(M+1,M+1,nDosSites,nBands*2);
    %    vkall=zeros(M+1,M+1,nDosSites,nBands);
        Ekall=zeros(M+1,M+1,nBands*2);
        else
        ukall=zeros(M+1,M+1,nDosSites,nBands);
        vkall=zeros(M+1,M+1,nDosSites,nBands);
        Ekall=zeros(M+1,M+1,nBands);
        end;
    end;
end;
% split casestring from directories
[~,cs1,cs2] = fileparts(casestring);
casestring=[cs1,cs2];
if ~ischar(Vimp)
    LDOSfileName0 = [casestring,'_Vimp_', num2str(Vimp),  '_N_', num2str(N)];
else
    % remove path from Vimp
    position=findstr(Vimp,filesep);
    if ~isempty(position)
        vimp=Vimp(position(length(position))+1:length(Vimp));
    else
        vimp=Vimp;
    end;
    LDOSfileName0 = [casestring,'_Vimp_', vimp,  '_N_', num2str(N)];
end;
% create some sub-directory to avoid many files in one directory
dirstring=[dirprefix,'data_',LDOSfileName0];
if division>0
    if ~(exist(dirstring,'dir'))
        mkdir(dirstring)
    end;
end;
% use external script to set up impurity Hamiltonian
HImpurity=get_Himp(Vimp,N,nOrbitals,sublattice);

% only one for loop
for index=startindex:endindex
    iKy= mod(index-1,M)+1;
    iKx= ceil(index/M);
    tic;
    ekukvk_file=[dirstring,'/','kx_',num2str(kx_ind(1,iKx)),'_',num2str(kx_ind(2,iKx)),'ky_',num2str(ky_ind(1,iKy)),'_',num2str(ky_ind(2,iKy)),'_g.mat'];
    ekukvk_fileGF=[dirstring,'/','kx_',num2str(kx_ind(1,iKx)),'_',num2str(kx_ind(2,iKx)),'ky_',num2str(ky_ind(1,iKy)),'_',num2str(ky_ind(2,iKy)),'_g_GF.mat'];
    calculate=true;
    if exist(ekukvk_fileGF, 'file') % no need to calculate this k-point
        calculate=false;
    end;
    if (exist(ekukvk_file,'file') && ~calcGreens) % no need to calculate this k-point if only ldos should be calculated
        calculate=false;
    end;
    if division==0
        calculate=true;
    end;
    if calculate
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
        [eigVector, eigValue] = eig(kSpaceHamiltonian);
        [eigValueK, sortingIndex] = sort(real(diag(eigValue)));
        if debug
            whos;
        end;
        clear eigValue
        eigVectorK = (eigVector(:,sortingIndex));
        if version >0
            Ek_vector=eigValueK;
        else
            Ek_vector=eigValueK((nBands + 1):end);
        end;
        clear eigValueK
        if ~calcGreens
            if version >0
                % only use this part for the calculation of the Greens
                % function
                uK = eigVectorK(siteIndices,:);
            else
            uK = eigVectorK(siteIndices,(nBands + 1):end);
            vK = eigVectorK(nBands + siteIndices,(nBands + 1):end);
            end
            clear eigVectorK
            % depending on the mode do different things
            if division==0
                % single calculation of full DOS
                Ek = repmat(Ek_vector, 1, nEnergyPoints) ;
                if version >0
                    greensKSpace(iKx, iKy, :, :) = ((abs(uK)).^2)*(1./(E - Ek + 1i*ita ));
                else
                    greensKSpace(iKx, iKy, :, :) = ((abs(uK)).^2)*(1./(E - Ek + 1i*ita )) + ...
                    ((abs(vK)).^2)*(1./(E + Ek + 1i*ita ));
                end
            else
                % save result in one single file
                % put k_vector in filename to avoid double calculation ?
                if version >0
                    save(ekukvk_file,'uK','Ek_vector');
                else
                    save(ekukvk_file,'uK','vK','Ek_vector');
                end           
            end
        else
            if version >0
                uK = eigVectorK(1:nBands,:);
            else
                uK = eigVectorK(1:nBands,(nBands+1):end);
                vK = eigVectorK((nBands+1):end,(nBands+1):end);
            end;
            if division==0
                EnRep = repmat(Ek_vector',nBands,1);
                if version >0
                    latticeGreensK(iKx, iKy, :, :) = (uK./(E - EnRep + 1i*ita))*(uK');
                else
                    latticeGreensK(iKx, iKy, :, :) = (uK./(E - EnRep + 1i*ita))*(uK') + (vK./(E + EnRep + 1i*ita))*(vK');
                end
            else
                if version >0
                save(ekukvk_fileGF,'uK','Ek_vector');
                else
                save(ekukvk_fileGF,'uK','vK','Ek_vector');
                end
            end;

        end
    else
        disp([ekukvk_file,' already calculated, skipping.' ])
    end
    toc;
end
for en=1:le
if ~calcGreens
if ~tetra
    LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_ita_', num2str(ita),sqstring];
else
    LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_tetra_corr']
end;
else
    LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_ita_', num2str(ita),'_e_',num2str(E(en)),sqstring];
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
                if debug
                    whos;
                end;
            end;
            try
                    ekukvk_file=[dirstring,'/','kx_',num2str(kx_ind(1,iKx)),'_',num2str(kx_ind(2,iKx)),'ky_',num2str(ky_ind(1,iKy)),'_',num2str(ky_ind(2,iKy)),'_g.mat'];
                    ekukvk_fileGF=[dirstring,'/','kx_',num2str(kx_ind(1,iKx)),'_',num2str(kx_ind(2,iKx)),'ky_',num2str(ky_ind(1,iKy)),'_',num2str(ky_ind(2,iKy)),'_g_GF.mat'];
                    if calcGreens
                        load(ekukvk_fileGF);
                        if version >0
                        else
                        if ~(exist('vK','var'))
                            vK=0*uK;
                        end;
                        end
                    else
                        try
                            load(ekukvk_file);
                        catch exception
                            clear uK vK
                            load(ekukvk_fileGF);
                            uK = uK(siteIndices,:);
                            if version >0
                            if ~(exist('vK','var'))
                                vK=0*uK;
                            else
                                vK = vK(siteIndices,:);
                            end
                            else
                                vK = vK(siteIndices,:);
                            end
                        end
                        if version >0
                        else
                        if ~(exist('vK','var'))
                            vK=0*uK;
                        end;
                        end
                    end;
            catch exception
                % missing k-point (or wrong input as number of k-points)
                % First can happen if one job crashes; catch this by
                % calculating on the fly
                disp('Missing k-point, recalculating on the fly.')
                impurity_dos(inputfile, division, -index);
                if calcGreens
                    load(ekukvk_fileGF);
                else
                    load(ekukvk_file);
                end;

            end
            % careful: double code here, change both when doing any
            % modifications
            if ~tetra
                if ~calcGreens
                    Ek = repmat(Ek_vector, 1, nEnergyPoints) ;
                    if version>0
                     greensKSpace(iKx, iKy, :, :) = ((abs(uK)).^2)*(1./(E(en) - Ek + 1i*ita ));
                    else
                    greensKSpace(iKx, iKy, :, :) = ((abs(uK)).^2)*(1./(E(en) - Ek + 1i*ita )) + ...
                        ((abs(vK)).^2)*(1./(E + Ek + 1i*ita ));
                    end;
                else
                    EnRep = repmat(Ek_vector',nBands,1);
                    if version>0
                    latticeGreensK(iKx, iKy, :, :) = (uK./(E(en) - EnRep + 1i*ita))*(uK');
                    else
                    latticeGreensK(iKx, iKy, :, :) = (uK./(E(en) - EnRep + 1i*ita))*(uK') + (vK./(E(en) + EnRep + 1i*ita))*(vK');
                    end;
                end;
            else
                if version >0
                ukall(iKx,iKy,:,:)=uK;
                else
                 ukall(iKx,iKy,:,:)=uK;
                vkall(iKx,iKy,:,:)=vK;
                end;
                Ekall(iKx,iKy,:)=Ek_vector;
                % store the edges twice to construct set of triangles that cover the whole area
%                 if iKx==1
%                    ukall(M+1,iKy,:,:)=uK;
%                    vkall(M+1,iKy,:,:)=vK;
%                    Ekall(M+1,iKy,:,:)=Ek_vector;
%                 end;
%                 if iKy==1
%                    ukall(iKx,M+1,:,:)=uK;
%                    vkall(iKx,M+1,:,:)=vK;
%                    Ekall(iKx,M+1,:,:)=Ek_vector;
%                 end;
%                 if (iKx==1) && (iKy==1)
%                    ukall(M+1,M+1,:,:)=uK;
%                    vkall(M+1,M+1,:,:)=vK;
%                    Ekall(M+1,M+1,:,:)=Ek_vector;
%                 end;
            end;
        end;
        disp('... done.');
    end
    % do the calculation of dos
    if ~calcGreens
        greensRealSpace = zeros(nDosSites, nEnergyPoints);
        disp('calculating GF in real space...');
        if ~tetra
        for iSite = 1: nDosSites
            disp(['Done ',num2str(iSite),' of ', num2str(nDosSites), 'nDosSites']);
            %if ~tetra
                for iEnergyPoint = 1:nEnergyPoints
                    if singular_quad
                        greensRealSpace(iSite, iEnergyPoint) = (1/(2*pi))^2*delKx*delKy*singular_double_quad(1./squeeze(greensKSpace(:, :, iSite, iEnergyPoint)));
                    else
                        greensRealSpace(iSite, iEnergyPoint) = (1/(2*pi))^2*delKx*delKy*sum(sum(squeeze(greensKSpace(:, :, iSite, iEnergyPoint))));
                    end;
                end
        end
            %else
        else
                disp('...using 2D version of Tetrahedron method');
                % set up a k-mesh that is suitable to cover the whole
                % Brillouinzone with triangles
                mesh1=[0.5:1:(M+0.5)]*2*pi/M;
                % to do: kx,ky can be only a vector to simplify indexing
                [kx,ky] = meshgrid(mesh1, mesh1);
                    % to do: vectorize the code!
                    for iband=1:nBands
                        disp(['Band ',num2str(iband),' of ',num2str(nBands)]);
                        if version >0
                            E=Ekall(:,:,iband);
                            a=ukall(:,:,:,iband).*conj(ukall(:,:,:,iband));
                            greensRealSpace(:, :) = greensRealSpace(:, :) + f(E,a,kx,ky,energy);
                            E=Ekall(:,:,iband+nBands);
                            a=ukall(:,:,:,iband+nBands).*conj(ukall(:,:,:,iband+nBands));
                            greensRealSpace(:, :) = greensRealSpace(:, :) + f(E,a,kx,ky,energy);                           
                        else
                        E=Ekall(:,:,iband);
                        a=ukall(:,:,:,iband).*conj(ukall(:,:,:,iband));
                        greensRealSpace(:, :) = greensRealSpace(:, :) + f(E,a,kx,ky,energy);
                        a=vkall(:,:,:,iband).*conj(vkall(:,:,:,iband));
                        greensRealSpace(:, :) = greensRealSpace(:, :) + f(E,a,kx,ky,-energy);
                        end;
                    end;
        end;
    else
        latticeGreens = zeros(nBands, nBands);
        for i = 1:nBands
            if debug
                whos;
            end;
            disp(['Integrating Bands (',num2str(i),' ,:) of ', num2str(nBands),'.']);
            for j = 1:nBands
                if singular_quad
                    latticeGreens(i, j) = (1/(2*pi))^2*delKx*delKy*singular_double_quad(1./squeeze(latticeGreensK(:,:,i,j)));
                else
                    latticeGreens(i, j) = (1/(2*pi))^2*delKx*delKy*sum(sum(squeeze(latticeGreensK(:,:,i,j))));
                end
            end
        end
    end
        %else
%     if tetra % obsolete (already above!)
%             disp('...using 2D version of Tetrahedron method');
%             % set up a k-mesh that is suitable to cover the whole
%             % Brillouinzone with triangles
%             mesh1=[0:1:(M)]*2*pi/M;
%             % to do: kx,ky can be only a vector to simplify indexing
%             [kx,ky] = meshgrid(mesh1, mesh1);
%                 % to do: vectorize the code!
%                 for iband=1:nBands
%                     disp(['Band ',num2str(iband),' of ',num2str(nBands)]);
%                     E=Ekall(:,:,iband);
%                     a=ukall(:,:,:,iband).*conj(ukall(:,:,:,iband));
%                     greensRealSpace(:, :) = greensRealSpace(:, :) + f(E,a,kx,ky,energy);
%                     a=vkall(:,:,:,iband).*conj(vkall(:,:,:,iband));
%                     greensRealSpace(:, :) = greensRealSpace(:, :) + f(E,a,kx,ky,-energy);
%                 end;
%     end;
      %  end;
   % end
    disp('Writing out LDOS ...');
    if isempty(dirprefix)
        outputfilename=[LDOSfileName];
    else
        outputfilename=[dirprefix,LDOSfileName];
    end;

    if ~calcGreens
        if tetra
            ldos=greensRealSpace;
        else
            ldos = (-(1/pi))*imag(greensRealSpace);
        end;
        if isempty(nDosSitesfile)
            orbitalLDOSFarAway = ldos(1:efforb,:);
            totalLDOSFarAway = sum(orbitalLDOSFarAway,1); %#ok<NASGU>
            orbitalLDOSImp = ldos(efforb+1:2*efforb,:);
            totalLDOSImp = sum(orbitalLDOSImp,1);%#ok<NASGU>
            orbitalLDOSImpNN = ldos(2*efforb+1:3*efforb,:);
            totalLDOSImpNN = sum(orbitalLDOSImpNN,1);%#ok<NASGU>
            orbitalLDOSImpNNN = ldos(3*efforb+1:4*efforb,:);
            totalLDOSImpNNN = sum(orbitalLDOSImpNNN,1);%#ok<NASGU>
            % to be done: change filename to general string
            save(outputfilename, 'energy', 'orbitalLDOSFarAway', 'orbitalLDOSImp', 'orbitalLDOSImpNN', 'orbitalLDOSImpNNN');
        else
            % new output format
            for n=1:nDosSites/efforb
                orbitalLDOS((n-1)*efforb+(1:efforb),:)=ldos((n-1)*efforb+(1:efforb),:);
            end;
             save(outputfilename, 'energy', 'orbitalLDOS', 'LDOSsites','efforb');           
        end
    else
        % some arkward workaround for multiple energies
        E1=E;
        E=E1(en);
        save(outputfilename,'latticeGreens','N','nOrbitals','E','sublattice');
        E=E1;
    end;
end
end;
r=1;
