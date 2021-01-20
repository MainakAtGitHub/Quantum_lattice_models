function r=BdG_impurity_v3(inputfile,mode,saveEnTot);%,normal_metal)%,dislocation_length)%,pos_file)

% Modified impurity BdG code to include
% 1. Convergence check parameter as 
%       (a) nCheck = abs(nCal - n)/n;
%       (c) deltaCheck = norm(deltaCal - delta)/norm(delta);
% 2. Randomize mixing parameter
%       First find beta for converging solution. Now choose a range close
%       to this beta, say [beta1 beta2] and for each iteration take new
%       beta to be beta = rand
global fullgamma cutek dress saveEigVec BdGfileName

if isempty(saveEigVec)
    saveEigVec = false;
end

% Mainak
% normal_metal = true;
% normal_metal = false;
% Mainak

% if nargin<5
%     dislocation_length=false;
% end
% Mainak
% if nargin< 4
%     normal_metal = false;
% end
% Mainak
if ~exist('normal_metal','var')
    normal_metal=false;
end
% MAINAK
if nargin < 3
    saveEnTot = false;
end
% MAINAK
if nargin < 2
    mode=0;
else
    if isdeployed
        mode=str2double(mode);
    end
end
% set some default value (kept for backward compatibility)
% fast summation with energies close to 0
cutek=NaN;
if nargin <1
    %default Parameters
    N = 9;%input('Enter N   ');
    Vimp = .4;%input('enter impurity potential    ');
    alpha = .25;%input('enter alpha    '); % self-consistency parametser
    beta1 = .6;%input('enter beta1    ');
    beta2 = .9;%input('enter beta1    ');
    deltaTol = 1e-4;%input('enter tolerance for gap convergence    ');
    nTol = 1e-3;%input('enter tolerance for electron density convergence    ');
    maxLoop = 60;%input('enter maxloop     '); % max no of iterations for self consistency
    nOrbitals = 10;
    n0 = 1.2*nOrbitals; % no. of valence electrons per unit cell
    kT = .01;
    TB_file='TB_hamiltonian_FeSe_2D.mat';
    Gamma_file='Gamma_FeSe_Toms_BS_6Dec13_cut_2.mat';
    BdGfileName = ['BdG_Impurity_FeSe', '_N_', num2str(N),'_Vimp_', num2str(Vimp)];
    %  casestring='LDOS_FeSe_Milan_Gamma';
    input_fileName = ['BdG_homogeneous_FeSe_Toms_BS_6Dec13', '_N_', num2str(N),'_GammaCut_',num2str(2),'.mat'];
    input_fileName = ['BdG_Impurity_FeSe_Toms_BS_6Dec13_N_9_Vimp_0.4']
    BdGfileName = ['BdG_Impurity_FeSe_Toms_BS_6Dec13', '_N_', num2str(N),'_Vimp_', num2str(Vimp)];
else
    % otherwise read inputfile
    try
        % old input format with mat-file (backward compatibility)
        load(inputfile);
    catch %err
        % new text-based input format (reads in text file and sets
        % variables accordingly)
        read_input_file=inputfile;
        read_input;
        read_input_file
    end
end


% Mainak

if ~exist('bilayer_int','var')
    bilayer_int=false;
end
if ~exist('int_soc','var')
    int_soc=false;
end

% Mainak


% input files: tight-binding model
load(TB_file)
% any reason for these double variables?
% latticeVectors = latticeVector;
% pairing interaction
load(Gamma_file,'-mat')
% default: same input as output filename, can be removed later
%if ~(exist('input_fileName','var'))
%    input_fileName=BdGfileName;
%end;
%%%%%Mainak
if exist('ref_grid_hopping_file','var')
    t_r_ref=load(ref_grid_hopping_file);
end
%%%%%Mainak
%%%% Mainak
% if nargin>5

if ~exist('pos_file','var')
    dislocation_length=0
else
    load(pos_file,'-mat');
    dislocation_length=N^2-size(r,1);
end
% end
%%%%Mainak
% load mean fields from previous iteration or seed
% fix for existing function mu
mu=0;
load(BdGfileName,'-mat');
[filepath,name,ext]=fileparts(BdGfileName);
if ~exist('sublattice','var')
    % sublattice= {-1,0,1} to define whether there are two sites per
    % elementary cell and which site is first
    sublattice=1
end
% a switch to produce more output
if ~exist('debug','var')
    debug=false;
end
% condition to stop cycle if gap converges to zero
if ~exist('deltaTol','var')
    deltaTol=eps;
end
% switch to activate memory management actions (clear, sparse matrix
% arrays)
if ~exist('memorymanagement','var')
    memorymanagement=false;
end
% magnetic calculation ?
if ~exist('magnetic','var')
    magnetic=false;
end
% supercell calculation
if ~exist('super','var')
    super=false;
end
% supercell size for BdG
if (~exist('M_super','var') && (super== true))
    M_super=M;
end
% spin polarized calculation (including a Zeeman magnetic field)?
if ~exist('spinpolarized','var')
    spinpolarized=false;
end

if ~exist('spinfullnormal','var')
    spinfullnormal=false;
end

% BdG matrix blocks (to be done: make it work for non-square system sizes)
%%%%% Mainak temporary changing size due to dislocation
if true %dislocation_length>0
    nBands = (N^2-dislocation_length)*nOrbitals;
else
    nBands = N^2*nOrbitals;
end
%%%%% Mainak

if ~super
    % kinetic energy
    % Mainak
    if true %dislocation_length>0
        if ~exist('eff_pot','var')
            eff_pot=1.5;
        end
        if ~exist('further_N_cut_off','var')
            further_N_cut_off=1.45;
        end
        if ~exist('NN_cut_off','var')
            NN_cut_off=1.3;
        end
        if exist('ref_grid_hopping_file','var') 
%             if nOrbitals==2
                [H0,H_soc1,H_soc2] = hoppings(r,N,true,t_r_ref,eff_pot,further_N_cut_off,NN_cut_off,filepath,nOrbitals,bilayer_int,int_soc);
%             else
%                 [H0,H_soc1,H_soc2] = hoppings(r,N,true,t_r_ref,eff_pot,further_N_cut_off,NN_cut_off,filepath,nOrbitals,false);
%             end
        else 
            [H0,H_soc1,H_soc2] = hoppings(r,N,true,[],eff_pot,further_N_cut_off,NN_cut_off);
        end
%         H0 = plot_line_disloc_hor((N-1)/2,dislocation_length);
    else
        % Mainak
        H0 = lattice_translation(N, TBparameters, latticeVector);
    end
else
    [H0, superLatticeVectors] = supercell_hoppings(N, TBparameters, latticeVector);
end

% ugly global variable to treat full orbital dependent pairing interaction
if exist('Gamma','var')
    % Mainak%%%% Setting SCInteractionMatrix=0 by hand now, can probably be
    % generated by hoppings() function later
    if normal_metal
        SCInteractionMatrix = zeros(size(H0));
    else
        % Mainak
        SCInteractionMatrix = lattice_translation(N, Gamma, latticeVectorsSC);
    end
    fullgamma=false;
else
    fullgamma=true;
    SCInteractionMatrix=Gammafull;
end

% fix the position of the impurity unit cell
% to be done: allow for different impurity positions; allow for multiple
% impurities
impCell = [ceil(N/2) ceil(N/2)];
% % allow for general impurity potentials
% if ~ischar(Vimp)
%     [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, impCell, impCell);
% else
%     load(Vimp,'-mat')
%     % now we have a set of impurity matrices imp_matr
%     % together with some lattice vectors imp_vec
%     % here we set up Himp directly
%     numimp=size(imp_vec,2);
%     for n=1:numimp
%         cellvector=impCell+imp_vec(:,n);
%         cellvector(1)=mod(cellvector(1)+ceil(N/2),N);
%         cellvector(2)=mod(cellvector(2)+ceil(N/2),N);
%         [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, impCell, cellvector);
%         Himp(iRange, jRange)=imp_matr(:,:,n);
%     end;
% end;
impNNCell = impCell + [0 1];
[iImpNNRange, jImpNNRange] = find_lattice_translation_index(N, nOrbitals, impNNCell, impCell);

switch sublattice
    case 1
%         % default with first iron in right upper corner (FeSe)
%         % allow for different potentials
%         if ~ischar(Vimp)
%         if numel(Vimp)==1
%             impPotential = Vimp*eye(nOrbitals/2, nOrbitals/2);
%         else
%             impPotential = diag(Vimp);
%         end;
%         % to do: setup correct impurity potential
%         Himp(iRange, jRange) = [impPotential zeros(nOrbitals/2); zeros(nOrbitals/2) zeros(nOrbitals/2)];
%         else
%         end;
        
        % Indices of sites NN and NNN to impurity (for status output only)
        iNNsiteRange = iImpNNRange(1:nOrbitals/2);
        jNNsiteRange = jImpNNRange(1:nOrbitals/2);
        iNNNsiteRange = iImpNNRange((1+nOrbitals/2):nOrbitals);
        jNNNsiteRange = jNNsiteRange;
    case -1
        % to be implemented (not tested), reversed relative position of
        % Fe(1) and Fe(2)
        iNNsiteRange = iImpNNRange(1:nOrbitals/2);
        jNNsiteRange = jImpNNRange(1:nOrbitals/2);
        iNNNsiteRange = iImpNNRange((1+nOrbitals/2):nOrbitals);
        jNNNsiteRange = jNNsiteRange;
    case 0
%         % no sublattice (5 orbital or 1 band model)
%         if ~ischar(Vimp)
%         if numel(Vimp)==1
%             Himp(iRange, jRange) = Vimp*eye(nOrbitals, nOrbitals);
%         else
%             Himp(iRange, jRange) = diag(Vimp);
%         end;   
%          else
%         end;
        iNNsiteRange = iImpNNRange(1:nOrbitals);
        jNNsiteRange = jImpNNRange(1:nOrbitals);
        impNNNCell = impCell + [1 1];
        [iImpNNNRange, jImpNNNRange] = find_lattice_translation_index(N, nOrbitals, impNNNCell, impCell);
        iNNNsiteRange = iImpNNNRange(1:nOrbitals);
        jNNNsiteRange = jImpNNNRange(1:nOrbitals);
end
%Self consistency iteration
% initial guess for the density of up-electrons (if not set in the input
% already)
if ~(exist('nUp','var'))
    fillup=0.5*n0/nOrbitals;
    nUp = fillup*ones(nBands,1);
end
% same for down electrons (set to half filling as well)
if ~(exist('nDown','var'))
    filldown=0.5*n0/nOrbitals;
    nDown = filldown*ones(nBands,1);    
end

% Mainak
if ~(exist('uS','var'))
    uS = zeros(size(nUp));
end
if ~(exist('vS','var'))
    vS = zeros(size(nDown));
end
% Mainak

% Mainak
if ~(exist('U','var'))
    U = 0;
end
if ~(exist('U_pr','var'))
    U_pr = 0;
end
if ~(exist('J','var'))
    J = 0;
end
% Mainak


% % Mainak
% if spinpolarized
%     n_randdd = 0.01*rand(length(nUp),1);
%     nUp = nUp+n_randdd;
% end
% % Mainak
% 
% % Mainak
% if spinpolarized
%     nDown = nDown-n_randdd;
% end
% % Mainak


%%%%%%%%%%%%%%%% Mainak
%         nUpdown = nDown;
%         nDowndown = nUp;
%%%%%%%%%%%%%%%% Mainak

% New stuff here: correlated electrons (including U, U', J, J')
if ~(exist('correlated','var'))
    correlated=false;
end
if correlated
    spinpolarized=true;
end

if spinpolarized
    if ~(exist('nUpdown','var'))
        nUpdown = nDown;
    end;
    if ~(exist('nDowndown','var'))
        nDowndown = nUp;
    end;
end
%%%%%%%%Mainak
if spinfullnormal
    if ~(exist('nAnoUpDown','var'))
        nAnoUpDown=zeros(size(nUp));
%         nAnoUpDown=0.01*(rand(size(nUp))+1i*rand(size(nUp)))/sqrt(2);
    end
    if ~(exist('nAnoDownUp','var'))
        nAnoDownUp=conj(nAnoUpDown);
    end

end
%%%%%Mainak
% setup of some "growing" variables

% not needed for long time, remove
%if ~(exist('nUpAcc','var'))
% nUpAcc = [];
%end;
%if ~(exist('nDownAcc','var'))
% nDownAcc = [];
%end;
if ~(exist('deltaMaxAcc','var'))
 deltaMaxAcc = [];
end
if ~(exist('deltaMinAcc','var'))
 deltaMinAcc = [];
end
if ~(exist('deltaDiffAcc','var'))
 deltaDiffAcc = [];
end
if ~(exist('muAcc','var'))
 muAcc = [];
end
if ~(exist('nAcc','var'))
 nAcc=[];
end
if ~(exist('nUpAcc','var'))
 nUpAcc=[];
end
% by default mix delta
if ~(exist('mixdelta','var'))
    mixdelta=true;
end
if ~(exist('hom','var'))
    hom=false;
end
if ~(exist('dress','var'))
    dress=[];
end


% write out a warning
if ~mixdelta
    disp('Warning: Not mixing delta, only converging nUp, nDown, mu.');
end

% setting of Hamiltonian
Himp=get_Himp(Vimp,N,nOrbitals,sublattice,false,BdGfileName,dislocation_length);
if ~super
    H = H0 + Himp;
else
    [HSuper, superLatticeVectors] = supercell_hoppings(N, TBparameters, latticeVector);
end


if magnetic
    % do a magnetic simulation with magnetic impurity
    % some default behavior: If Vimpdown not defined, use a potential
    % scatterer
    if ~exist('Vimpdown','var')
        Vimpdown=Vimp
    end
    % set the impurity Hamiltonian for down spins
    Himpdown=get_Himp(Vimpdown,N,nOrbitals,sublattice);
    % set the normal state Hamiltonian for down electrons:
    Hdown=H0+Himpdown;
    Hdown=0.5*(Hdown+Hdown'); % make it explicitely Hermitian (numerical inaccurancy)
    % save memory
    clear Himpdown
    % if spin-polarized, i.e. with Zeeman term, use two chemical potentials
    % for up and down electrons
    if spinpolarized
        if ~exist('mudown','var')
            mudown=mu
        end
    end
end

% avoid some numerical inaccurancy; for some reason lattice_translation as
% well as the general impurity potential (hoppings)
% gives back a non-hermitian matrix with sum(sum(abs(H-H'))) ~ 1e-13
if ~super
    H=0.5*(H+H');
    clear Himp;
    clear H0;
else
    %nSuperCells = size(superLatticeVectors,1);
    %for iUnitCell = 1:nSuperCells
    %    HSuper(:,:,iUnitCell)=0.5*(HSuper(:,:,iUnitCell)+HSuper(:,:,iUnitCell)');
    %end;
end
% single shot calculation: just diagonalize and exit
if mode==1
    maxLoop=1;
end

% BdG iterations: up to maxLoop

%Mainak
IntUp1 = zeros(size(nUp));
IntDown1 = zeros(size(nDown));
%Mainak

for i = 1:maxLoop
    % two cases for supercell calculations here, first the usual one
    if ~super
        % setup of the kinetic energy (including chemical potential)
        KE = H - mu*eye(nBands);
        
        if ~spinpolarized
            if ~magnetic
                [ nUpCal, nDownCal, deltaCal, En, TotKE ] = BdG_step( KE,delta, kT,nBands, SCInteractionMatrix);
            else
                KEdown = conj(Hdown - mu*eye(nBands));
                [ nUpCal, nDownCal, deltaCal, En, TotKE ] = BdG_step( KE,delta, kT,nBands, SCInteractionMatrix,-KEdown);
            end
        else
            mudown=mu;
            % put a magnetic field here
            if exist('field','var')
%               mu=mu+field;
                mudown=mu-field(end);
                KE=KE-field(end)*eye(nBands);
%               mu=mu+field;
            end
            if ~magnetic
                %%%%%%%%%%%%% Mainak
                if correlated
%                     U = 2.5;
                    KE = KE + U*diag(nDown);
                    for it1 = 1:length(nDown)
                        IntUp1(it1) = U_pr*(sum(nDown([nOrbitals*(ceil(it1/nOrbitals)-1)+1:nOrbitals*ceil(it1/nOrbitals)]))-nDown(it1)) + (U_pr-J)*(sum(nUp([nOrbitals*(ceil(it1/nOrbitals)-1)+1:nOrbitals*ceil(it1/nOrbitals)]))-nUp(it1));
                    end
                    KE = KE + diag(IntUp1);
                end
                %%%%%%%%%%%%% Mainak

                KEdown = H - mudown*eye(nBands);
                %%%%%%%%%%%%% Mainak
                if correlated
%                     U = 2.5;
                    KEdown = KEdown + U*diag(nUp);
                    for it2 = 1:length(nUp)
                        IntDown1(it2) = U_pr*(sum(nUp([nOrbitals*(ceil(it2/nOrbitals)-1)+1:nOrbitals*ceil(it2/nOrbitals)]))-nUp(it2)) + (U_pr-J)*(sum(nDown([nOrbitals*(ceil(it2/nOrbitals)-1)+1:nOrbitals*ceil(it2/nOrbitals)]))-nDown(it2));
                    end
                    KEdown = KEdown + diag(IntDown1);                    
                end
                %%%%%%%%%%%%% Mainak
                if spinfullnormal
                    H_off_up = -U*diag(nAnoUpDown);
                    H_off_down = -U*diag(nAnoDownUp); %%%??????????????COULD BE SETA AS conj(H_off_up)?
                    if exist('field','var')
                        if numel(field)>1                          
                            H_off_up=H_off_up-field(1)*eye(nBands);
                            H_off_down=H_off_down-field(1)*eye(nBands);

                            H_off_up=H_off_up-(-1i)*field(2)*eye(nBands);
                            H_off_down=H_off_down-(1i)*field(2)*eye(nBands);
                        end
                    end 
                    
                    if int_soc
                        H_off_up=H_off_up+H_soc1+(-1i)*H_soc2;
                        H_off_down=H_off_down+H_soc1+(1i)*H_soc2;
                    end
                end

                if ~spinfullnormal
                    [ nUpCal, nDownCal, deltaCal, En, TotKE ] = BdG_step( KE,delta, kT,nBands, SCInteractionMatrix,-KEdown);
                else
                    [ nUpCal, nDownCal, deltaCal, En, TotKE, nAnoUpDownCal, nAnoDownUpCal] = BdG_step1( KE,delta, kT,nBands, SCInteractionMatrix,KEdown,H_off_up,H_off_down);
                end
            
                % [ nUpCaldown, nDownCaldown, deltaCaldown ] = BdG_step( KEdown, conj(-delta'), kT, nBands, SCInteractionMatrix,-KE);
            else
                KEdown = Hdown - mudown*eye(nBands);
                [ nUpCal, nDownCal, deltaCal, En, TotKE ] = BdG_step( KE,delta, kT,nBands, SCInteractionMatrix,-KEdown);
                % [ nUpCaldown, nDownCaldown, deltaCaldown ] = BdG_step( KEdown, conj(-delta'), kT, nBands, SCInteractionMatrix,-KE);
            end
        end

    else
        if ~spinpolarized
            if ~magnetic
                % set the supercell parameters
                BZ=supercell_parameters(M_super);
                % 1 supercell only
                %BZ.k=[0 0;0 0 ; 0 0];
                %BZ.weight=[1; 1;1]/3;
                %BZ.mu=mu;
                % first version of supercell calculation (not well tested)
                [ nUpCal, nDownCal, deltaCal, En ] = BdG_step_super( HSuper,delta, kT,nBands, SCInteractionMatrix,BZ,Himp,mu);
            else
                disp('not implemented')
            end
        else
            disp('not implemented')
        end
    end
%     BdGMatrix = [KE -delta; -delta' -KE];
%     if memorymanagement
%         clear KE;
%     end;
%     [eVector, eValue] = eig(BdGMatrix);
%     % save some memory for following commands (here we need to save three
%     % full arrays such that we get in MB:
%     % 3*(2*N^2*nOrbitals)^2*8/1024/1024 (3.6G for N=25, 470M for N=15)
%     %clear BdGMatrix
%     [En, sortIndex] = sort(real(diag(eValue)));
%     % save some memory for following commands
%     %clear eValue
%     eVector = eVector(:,sortIndex);
%     fermi = 1./(1 + exp(En/kT));
%     nUpCal = (abs(eVector(1:nBands,:)).^2)*fermi;
%     nDownCal = (abs(eVector((nBands + 1):end,:)).^2)*(1 - fermi);
%     deltaCal = SCInteractionMatrix.*((eVector(1:nBands,:)*(((eVector((nBands + 1):end,:))').*repmat(fermi,1,nBands))));
    %clear eVector
    % convergence criterium: norm (as defined for vector)
    if mode==1
        dlmwrite([inputfile,'_fill'],(1/(numel(nUp)/nOrbitals))*sum(nUpCal + nDownCal),'precision',10);
        return;
    end
    deltaDiff = norm(deltaCal(:) - delta(:))/norm(delta(:));
    
    % Mainak
    nUpDiff = norm(nUpCal(:) - nUp(:))/norm(nUp(:));
    % Mainak
    
    nDiff = abs((1/(numel(nUp)/nOrbitals))*sum(nUpCal + nDownCal) - n0)/n0;
  %  if spinpolarized
   %     % to be checked
   %     tmp=-conj(deltaCaldown');
   %     deltaDiff(2)= norm(tmp(:) - delta(:))/norm(delta(:));
   %     nDiff(2) = abs((1/N^2)*sum(nUpCaldown + nDownCaldown) - n0)/n0;
   %     clear tmp;
   % end;
   
   % Mainak
   if saveEnTot
       %        uS = (abs(eVector(1:nBands,floor(nBands/2)).^2));
       %        vS = (abs(eVector((nBands + 1):end,floor(nBands/2))).^2);
       %        uS = (abs(eVector(1:nBands,:)).^2);
       %        vS = (abs(eVector((nBands + 1):end,:)).^2);
%        if ~spinfullnormal
       E_Hub = U/4*sum((nUp+nDown).^2 - (nUp-nDown).^2);%-(mag_x.^2)-mag_y.^2+%sum_x.^2+sum_y.^2;
%        else
%            E_Hub = U/4*sum((nUp+nDown).^2 - (nUp-nDown).^2 ...
%            - (nAnoUpDown+nAnoDownUp).^2 - (1i*(nAnoUpDown-nAnoDownUp)).^2);
%        end
       if spinfullnormal
           E_Hub_Ano = -U/4*sum((nAnoUpDown+nAnoDownUp).^2 + (1i*(nAnoUpDown-nAnoDownUp)).^2);
%            if exist('field','var')
%                if numel(field)>1
%                    E_field_x = sum(-field(1)*(nAnoUpDown+nAnoDownUp));
%                    E_field_y = sum(-field(2)*(1i)*(nAnoUpDown-nAnoDownUp));
%                end
%            end
       end

       E_Hub_12 = sum(IntUp1.*nUp) + sum(IntDown1.*nDown);          
       if Gamma(:,:,1) == 0
           E_Sup =0;
       else
           E_Sup = -1/Gamma(:,:,1)*sum(sum(abs(delta).^2));
       end
       if spinfullnormal
%        TotEn = 1/numel(nUp)*(TotKE + E_Hub + E_Hub_12 + E_Sup + E_Hub_Ano);
         TotEn = 1/numel(nUp)*(TotKE - E_Hub + E_Hub_12 + E_Sup - E_Hub_Ano) + mu*n0;
%            if exist('field','var')
%                if numel(field)>1
%                  TotEn = TotEn + 1/numel(nUp)*(E_field_x + E_field_y);  
%                end
%            end
       else 
       TotEn = 1/numel(nUp)*(TotKE + E_Hub + E_Hub_12 + E_Sup);
       end
%      save([BdGfileName,'eigTotEn'],'uS','vS','TotEn');
       save([BdGfileName,'eigTotEnNoVec'],'TotEn','-ascii');
       disp('Total Energy =');
       disp(TotEn);
   end
   % Mainak
   if ~normal_metal
       if ((sum(nDiff) < numel(nDiff)*nTol) && ((sum(deltaDiff) < numel(deltaDiff)*deltaTol)|| (norm(delta(:))/(numel(nUp)/nOrbitals) < deltaTol)))
           break % go out of loop if self-consistency is achieved
       end
   else
       if ((sum(nDiff) < numel(nDiff)*nTol))
           break % go out of loop if self-consistency is achieved
       end
   end
        
    % (Mainak) norm(delta(:))/N^2 < deltaTol condition to be revoked for normal metal? 
    % homogenize calculation to get faster convergence without impurity
    if (hom && (Vimp==0))
        delta=homogenize_delta(delta,latticeVectorsSC,nOrbitals,N);
    end
    % update
    beta =  beta1 + (beta2 - beta1).*rand(1); 
    nUp = beta*nUp + (1-beta)*nUpCal;
    nDown = beta*nDown + (1-beta)*nDownCal;
    if spinfullnormal
        nAnoUpDown=beta*nAnoUpDown + (1-beta)*nAnoUpDownCal;
        nAnoDownUp=beta*nAnoDownUp + (1-beta)*nAnoDownUpCal;
    end
    % Mainak
%     uS = beta*uS + (1-beta)*uCal;
%     uS = uCal;
%     vS = beta*vS + (1-beta)*vCal;
%     vS = vCal;
    % Mainak
    
  %  if spinpolarized
  %      % to be checked
  %      nUpdown = beta*nUpdown + (1-beta)*nUpCaldown;
  %      nDowndown = beta*nDowndown + (1-beta)*nDownCaldown;
  %  end;
    % new variable for input file: mixdelta to only converge nUp, nDown, mu
    % with fixing delta (makes only sense if the initial guess for delta is
    % already good).
    if mixdelta
      %  if ~spinpolarized
            delta = beta*delta + (1-beta)*deltaCal;
       % else
        %    delta = beta*delta + 0.5*(1-beta)*(deltaCal-deltaCaldown');
       % end
    end
%     nAvg = (1/N^2)*(sum(nUp + nDown));
    nAvg = (1/(numel(nUp)/nOrbitals))*(sum(nUp + nDown));
    % Mainak
%     nUpAvg = (1/N^2)*(sum(nUp));
    nUpAvg = (1/(numel(nUp)/nOrbitals))*(sum(nUp));
    % Mainak
    %BM convention
    %    nAvg = (1/N^2)*(sum(nUpCal + nDownCal));

   % if spinpolarized
   %     nAvg(2) = (1/N^2)*(sum(nUpdown + nDowndown));
   % end
    %if ~spinpolarized
     %   mu=mu - alpha*(nAvg - n0);
    %else
        mu=mu - alpha*(mean(nAvg) - n0);
        %mudown=mudown - alpha*(nAvg(2) - n0)
    %end;
    nAcc = [nAcc; mean(nAvg)];
    
    % Mainak
    nUpAcc = [nUpAcc; mean(nUpAvg)];
    % Mainak
    
    % fix phase of delta (mostly not necessary, but always gives the same
    % result, largest gap set to be positive
   % [~, index]=max(abs(delta(:)));
   % delta=delta*exp(-1i*angle(delta(index)));
    % second possible observables
    deltaMaxNN = max(max(abs(delta(iNNsiteRange, jNNsiteRange))));
    deltaMaxNNN = max(max(abs(delta(iNNNsiteRange, jNNNsiteRange))));
    %deltaMaxAcc = [deltaMaxAcc; deltamax];
    %deltaMinAcc = [deltaMinAcc; min(min(real(delta)))]; 
    deltaMaxAcc = [deltaMaxAcc; deltaMaxNN];
    deltaMinAcc = [deltaMinAcc; deltaMaxNNN];
    % include the lowest 50 eigenenergies into the file
    if exist('energiesAcc','var')
        energiesAcc=[energiesAcc; En(nBands+1:nBands+50)'];
    else
        energiesAcc=[En(nBands+1:nBands+50)'];
    end;
    muAcc = [muAcc; mu];
    deltaDiffAcc = [deltaDiffAcc; sum(deltaDiff)];
    disp([num2str(i),' ndiff= ',num2str(nDiff), ' deltaDiff= ',num2str( deltaDiff), ' deltaMaxNN= ',num2str(deltaMaxNN)]);
    if ~spinpolarized
        if size(delta,1)>11000
            save(BdGfileName,'nAcc','nUpAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','energiesAcc','muAcc','mu', 'nUp','nDown','-v7.3');
        else
            save(BdGfileName,'nAcc','nUpAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','energiesAcc','muAcc','mu', 'nUp','nDown');
        end
    else
        if size(delta,1)>11000
            if spinfullnormal
              save(BdGfileName,'nAcc','nUpAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','energiesAcc','muAcc','mu','mudown','nUp','nDown','nUpdown','nDowndown','nAnoUpDown','nAnoDownUp','-v7.3');
            else 
              save(BdGfileName,'nAcc','nUpAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','energiesAcc','muAcc','mu','mudown','nUp','nDown','nUpdown','nDowndown','-v7.3');                
            end
        else
            if spinfullnormal
              save(BdGfileName,'nAcc','nUpAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','energiesAcc','muAcc','mu','mudown','nUp','nDown','nUpdown','nDowndown','nAnoUpDown','nAnoDownUp');
            else
              save(BdGfileName,'nAcc','nUpAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','energiesAcc','muAcc','mu','mudown','nUp','nDown','nUpdown','nDowndown');
            end
        end
    end
end
if i < maxLoop
    disp('Converged')
    % save the converged result
%  save(output_fileName,'nAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','muAcc','mu', 'deltaTol', 'nTol');
else
    disp('***********Not converged**********')
end
%r=1; % generate an error to stop the program
% plot (not needed any more)
% figure;
% subplot(2,2,1); plot(nAcc); title('nAcc'); axis('square');
% subplot(2,2,2); plot(muAcc); title('mu'); axis('square');
% subplot(2,2,3); plot(deltaMaxAcc); title('deltaMaxNN'); axis('square');
% subplot(2,2,4); plot(deltaMinAcc); title('deltaMaxNNN'); axis('square');
% figure; plot(deltaDiffAcc); title('Norm deltaDiff'); axis('square');



        

