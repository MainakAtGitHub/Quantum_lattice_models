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
    catch err
        % new text-based input format
        read_input_file=inputfile;
        read_input;
        read_input_file
    end;
end;



% input files
load(TB_file)
% any reason for these double variables?
latticeVectors = latticeVector;
load(Gamma_file,'-mat')
% default: same input as output filename, can be removed later
if ~(exist('input_fileName','var'))
    input_fileName=BdGfileName;
end;
load(input_fileName,'-mat');
% not really necessary?
%deltaH = delta; 
%muH = mu;
%clear delta mu;

% BdG matrix blocks
nBands = N^2*nOrbitals;
H0 = lattice_translation(N, TBparameters, latticeVectors);
SCInteractionMatrix = lattice_translation(N, Gamma, latticeVectorsSC);
Himp = zeros(size(H0));
impCell = [ceil(N/2) ceil(N/2)];
[iRange, jRange] = find_lattice_translation_index(N, nOrbitals, impCell, impCell);
% allow for different potentials
if numel(Vimp)==1
    impPotential = Vimp*eye(nOrbitals/2, nOrbitals/2);
else
    impPotential = diag(Vimp);
end;
% debug
Himp(iRange, jRange) = [impPotential zeros(nOrbitals/2); zeros(nOrbitals/2) zeros(nOrbitals/2)];

% Indices of sites NN and NNN to impurity
impNNCell = impCell + [0 1];
[iImpNNRange, jImpNNRange] = find_lattice_translation_index(N, nOrbitals, impNNCell, impCell);
iNNsiteRange = iImpNNRange(1:nOrbitals/2);
jNNsiteRange = jImpNNRange(1:nOrbitals/2);
iNNNsiteRange = iImpNNRange((1+nOrbitals/2):nOrbitals);
jNNNsiteRange = jNNsiteRange;

%Self consistency iteration

% % initial guess
if ~(exist('nUp','var'))
    nUp = .6*ones(nBands,1);
end;
if ~(exist('nDown','var'))
    nDown = .6*ones(nBands,1);
end;
% setup of some "growing" variables
if ~(exist('nUpAcc','var'))
 nUpAcc = [];
end;
if ~(exist('nDownAcc','var'))
 nDownAcc = [];
end;
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
% writ out a Warning
if ~mixdelta
    disp('Warning: Not mixing delta, only converging nUp, nDown, mu.');
end;

% setting of Hamiltonian
H = H0 + Himp;
% BdG iterations
for i = 1:maxLoop
    KE = H - mu*eye(nBands);
    BdGMatrix = [KE -delta; -delta' -KE];
    [eVector eValue] = eig(BdGMatrix);
    % save some memory for following commands (here we need to save three
    % full arrays such that we get in MB:
    % 3*(2*N^2*nOrbitals)^2*8/1024/1024 (3.6G for N=25, 470M for N=15)
    clear BdGMatrix
    [En,sortIndex] = sort(real(diag(eValue)));
    % save some memory for following commands
    clear eValue
    eVector = eVector(:,sortIndex);
    fermi = 1./(1 + exp(En/kT));
    nUpCal = (abs(eVector(1:nBands,:)).^2)*fermi;
    nDownCal = (abs(eVector((nBands + 1):end,:)).^2)*(1 - fermi);
    deltaCal = SCInteractionMatrix.*((eVector(1:nBands,:)*(((eVector((nBands + 1):end,:))').*repmat(fermi,1,nBands))));
    clear eVector
    % convergence criterium: norm (as defined for vector)
    deltaDiff = norm(deltaCal(:) - delta(:))/norm(delta(:));
    nDiff = abs((1/N^2)*sum(nUpCal + nDownCal) - n0)/n0;
    if (nDiff < nTol) && (deltaDiff < deltaTol)
       break % go out of loop if self-consistency is achieved
    end
    % update
    beta =  beta1 + (beta2 - beta1).*rand(1); 
    nUp = beta*nUp + (1-beta)*nUpCal;
    nDown = beta*nDown + (1-beta)*nDownCal;
    % new variable for input file: mixdelta to only converge nUp, nDown, mu
    % with fixing delta (makes only sense if the initial guess for delta is
    % already good).
    if mixdelta
        delta = beta*delta + (1-beta)*deltaCal;
    end;
    nAvg = (1/N^2)*(sum(nUp + nDown));
    mu = mu - alpha*(nAvg - n0);
    nAcc = [nAcc; nAvg];   
    % fix phase of delta (mostly not necessary, but always gives the same
    % result, largest gap set to be positive
    [deltamax,index]=max(abs(delta(:)));
    delta=delta*exp(-1i*angle(delta(index)));
    % second possible observables
    deltaMaxNN = max(max(abs(delta(iNNsiteRange, jNNsiteRange))));
    deltaMaxNNN = max(max(abs(delta(iNNNsiteRange, jNNNsiteRange))));
    %deltaMaxAcc = [deltaMaxAcc; deltamax];
    %deltaMinAcc = [deltaMinAcc; min(min(real(delta)))]; 
    deltaMaxAcc = [deltaMaxAcc; deltaMaxNN];
    deltaMinAcc = [deltaMinAcc; deltaMaxNNN]; 
    muAcc = [muAcc; mu];
    deltaDiffAcc = [deltaDiffAcc; deltaDiff];
    disp([num2str(i),' ndiff= ',num2str(nDiff), ' deltaDiff= ',num2str( deltaDiff)]);  
    save(BdGfileName,'nAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','muAcc','mu', 'deltaTol', 'nTol','nUp','nDown');
end
if i < maxLoop
    disp('Converged')
    % save the converged result
%  save(output_fileName,'nAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','muAcc','mu', 'deltaTol', 'nTol');
else
    disp('***********Not Converged**********')
end

% plot
figure;
subplot(2,2,1); plot(nAcc); title('nAcc'); axis('square');
subplot(2,2,2); plot(muAcc); title('mu'); axis('square');
subplot(2,2,3); plot(deltaMaxAcc); title('deltaMax'); axis('square');
subplot(2,2,4); plot(deltaMinAcc); title('deltaMin'); axis('square');
figure; plot(deltaDiffAcc); title('Norm deltaDiff'); axis('square');



        

