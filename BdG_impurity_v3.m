function r=BdG_impurity_v3(inputfile)

% Modified impurity BdG code to include
% 1. Convergence check parameter as 
%       (a) nCheck = abs(nCal - n)/n;
%       (c) deltaCheck = norm(deltaCal - delta)/norm(delta);
% 2. Randomize mixing parameter
%       First find beta for converging solution. Now choose a range close
%       to this beta, say [beta1 beta2] and for each iteration take new
%       beta to be beta = rand

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
        % old input format with mat-file
        load(inputfile);
    catch %err
        % new text-based input format
        read_input_file=inputfile;
        read_input;
        read_input_file
    end;
end;



% input files
load(TB_file)
% any reason for these double variables?
% latticeVectors = latticeVector;
load(Gamma_file,'-mat')
% default: same input as output filename, can be removed later
%if ~(exist('input_fileName','var'))
%    input_fileName=BdGfileName;
%end;
load(BdGfileName,'-mat');
if ~exist('sublattice','var')
    % sublattice= {-1,0,1} to define whether there are two sites per
    % elementary cell and which site is first
    sublattice=1
end;
% a switch to produce more output
if ~exist('debug','var')
    debug=false;
end;
% switch to activate memory management actions (clear, sparse, 
if ~exist('memorymanagement','var')
    memorymanagement=false;
end;
% magnetic calculation
if ~exist('magnetic','var')
    magnetic=false;
end;
% spin polarized calculation
if ~exist('spinpolarized','var')
    spinpolarized=false;
end;
% BdG matrix blocks
nBands = N^2*nOrbitals;
% kinetic energy
H0 = lattice_translation(N, TBparameters, latticeVector);
SCInteractionMatrix = lattice_translation(N, Gamma, latticeVectorsSC);
%Himp = zeros(size(H0));
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
        
        % Indices of sites NN and NNN to impurity
        iNNsiteRange = iImpNNRange(1:nOrbitals/2);
        jNNsiteRange = jImpNNRange(1:nOrbitals/2);
        iNNNsiteRange = iImpNNRange((1+nOrbitals/2):nOrbitals);
        jNNNsiteRange = jNNsiteRange;
    case -1
        % to be implemented (not tested)
        iNNsiteRange = iImpNNRange(1:nOrbitals/2);
        jNNsiteRange = jImpNNRange(1:nOrbitals/2);
        iNNNsiteRange = iImpNNRange((1+nOrbitals/2):nOrbitals);
        jNNNsiteRange = jNNsiteRange;
    case 0
%         % no sublattice (5 orbital)
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
end;
%Self consistency iteration
% % initial guess
if ~(exist('nUp','var'))
    fillup=0.5*n0/nOrbitals;
    nUp = fillup*ones(nBands,1);
end;
if ~(exist('nDown','var'))
    filldown=0.5*n0/nOrbitals;
    nDown = filldown*ones(nBands,1);
end;
if spinpolarized
    if ~(exist('nUpdown','var'))
        nUpdown = nDown;
    end;
    if ~(exist('nDowndown','var'))
        nDowndown = nUp;
    end;
end
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
end;
if ~(exist('deltaMinAcc','var'))
 deltaMinAcc = [];
end;
if ~(exist('deltaDiffAcc','var'))
 deltaDiffAcc = [];
end;
if ~(exist('muAcc','var'))
 muAcc = [];
end;
if ~(exist('nAcc','var'))
 nAcc=[];
end;
% by default mix delta
if ~(exist('mixdelta','var'))
    mixdelta=true;
end;
if ~(exist('hom','var'))
    hom=false;
end;
% write out a warning
if ~mixdelta
    disp('Warning: Not mixing delta, only converging nUp, nDown, mu.');
end;

% setting of Hamiltonian
Himp=get_Himp(Vimp,N,nOrbitals,sublattice);
H = H0 + Himp;
clear Himp;
if magnetic
    % do a non-magnetic simulation
    if ~exist('Vimpdown','var')
        Vimpdown=Vimp
    end;
    Himpdown=get_Himp(Vimpdown,N,nOrbitals,sublattice);
    Hdown=H0+Himpdown;
    clear Himpdown
    if spinpolarized
        if ~exist('mudown','var')
            mudown=mu
        end;
    end;
end;

% BdG iterations
for i = 1:maxLoop
    KE = H - mu*eye(nBands);
    if ~spinpolarized
        if ~magnetic
            [ nUpCal, nDownCal, deltaCal ] = BdG_step( KE,delta, kT,nBands, SCInteractionMatrix);
        else
            KEdown = Hdown - mu*eye(nBands);
            [ nUpCal, nDownCal, deltaCal ] = BdG_step( KE,delta, kT,nBands, SCInteractionMatrix,-KEdown);
        end;
    else
        mudown=mu;
        % put a magnetic field here
        if exist('field','var')
            mu=mu+field;
            mudown=mu-field;
        end;    
        if ~magnetic
            KEdown = H - mudown*eye(nBands);
            [ nUpCal, nDownCal, deltaCal ] = BdG_step( KE,delta, kT,nBands, SCInteractionMatrix,-KEdown);
           % [ nUpCaldown, nDownCaldown, deltaCaldown ] = BdG_step( KEdown, conj(-delta'), kT, nBands, SCInteractionMatrix,-KE);
        else
            KEdown = Hdown - mudown*eye(nBands);
            [ nUpCal, nDownCal, deltaCal ] = BdG_step( KE,delta, kT,nBands, SCInteractionMatrix,-KEdown);
           % [ nUpCaldown, nDownCaldown, deltaCaldown ] = BdG_step( KEdown, conj(-delta'), kT, nBands, SCInteractionMatrix,-KE);
        end;
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
    deltaDiff = norm(deltaCal(:) - delta(:))/norm(delta(:));
    nDiff = abs((1/N^2)*sum(nUpCal + nDownCal) - n0)/n0;
  %  if spinpolarized
   %     % to be checked
   %     tmp=-conj(deltaCaldown');
   %     deltaDiff(2)= norm(tmp(:) - delta(:))/norm(delta(:));
   %     nDiff(2) = abs((1/N^2)*sum(nUpCaldown + nDownCaldown) - n0)/n0;
   %     clear tmp;
   % end;
    if (sum(nDiff) < numel(nDiff)*nTol) && (sum(deltaDiff) < numel(deltaDiff)*deltaTol)
       break % go out of loop if self-consistency is achieved
    end
    % homogenize calculation to get faster convergence without impurity
    if (hom && (Vimp==0))
        delta=homogenize_delta(delta,latticeVectorsSC,nOrbitals,N);
    end
    % update
    beta =  beta1 + (beta2 - beta1).*rand(1); 
    nUp = beta*nUp + (1-beta)*nUpCal;
    nDown = beta*nDown + (1-beta)*nDownCal;
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
    end;
    nAvg = (1/N^2)*(sum(nUp + nDown));
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
    % fix phase of delta (mostly not necessary, but always gives the same
    % result, largest gap set to be positive
    [~, index]=max(abs(delta(:)));
    delta=delta*exp(-1i*angle(delta(index)));
    % second possible observables
    deltaMaxNN = max(max(abs(delta(iNNsiteRange, jNNsiteRange))));
    deltaMaxNNN = max(max(abs(delta(iNNNsiteRange, jNNNsiteRange))));
    %deltaMaxAcc = [deltaMaxAcc; deltamax];
    %deltaMinAcc = [deltaMinAcc; min(min(real(delta)))]; 
    deltaMaxAcc = [deltaMaxAcc; deltaMaxNN];
    deltaMinAcc = [deltaMinAcc; deltaMaxNNN];
    muAcc = [muAcc; mu];
    deltaDiffAcc = [deltaDiffAcc; sum(deltaDiff)];
    disp([num2str(i),' ndiff= ',num2str(nDiff), ' deltaDiff= ',num2str( deltaDiff), ' deltaMaxNN= ',num2str(deltaMaxNN)]);
    if ~spinpolarized
        save(BdGfileName,'nAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','muAcc','mu', 'nUp','nDown');
    else
        save(BdGfileName,'nAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','muAcc','mu','mudown','nUp','nDown','nUpdown','nDowndown');
    end
end
if i < maxLoop
    disp('Converged')
    % save the converged result
%  save(output_fileName,'nAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','muAcc','mu', 'deltaTol', 'nTol');
else
    disp('***********Not converged**********')
end
r=1;
% plot (not needed any more)
% figure;
% subplot(2,2,1); plot(nAcc); title('nAcc'); axis('square');
% subplot(2,2,2); plot(muAcc); title('mu'); axis('square');
% subplot(2,2,3); plot(deltaMaxAcc); title('deltaMaxNN'); axis('square');
% subplot(2,2,4); plot(deltaMinAcc); title('deltaMaxNNN'); axis('square');
% figure; plot(deltaDiffAcc); title('Norm deltaDiff'); axis('square');



        

