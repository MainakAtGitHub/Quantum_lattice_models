
% Parameters
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


% input files
load TB_hamiltonian_FeSe_2D.mat
latticeVectors = latticeVector;
load Gamma_FeSe_Toms_BS_6Dec13_cut_2
fileName = ['BdG_homogeneous_FeSe_Toms_BS_6Dec13', '_N_', num2str(N),'_GammaCut_',num2str(2),'.mat'];
load(fileName);
deltaH = delta; 
muH = mu;
clear delta mu;


% BdG matrix blocks
nBands = N^2*nOrbitals;
H0 = lattice_translation(N, TBparameters, latticeVectors);
SCInteractionMatrix = lattice_translation(N, Gamma, latticeVectorsSC);
Himp = zeros(size(H0));
impCell = [ceil(N/2) ceil(N/2)];
[iRange, jRange] = find_lattice_translation_index(N, nOrbitals, impCell, impCell);
impPotential = Vimp*eye(nOrbitals/2, nOrbitals/2);
Himp(iRange, jRange) = [impPotential zeros(nOrbitals/2); zeros(nOrbitals/2) zeros(nOrbitals/2)];


%Self consistency iteration

% initial guess
nUp = .6*ones(nBands,1);
nDown = .6*ones(nBands,1);
mu = muH;
delta = deltaH;
nUpAcc = [];
nDownAcc = [];
deltaMaxAcc = [];
deltaMinAcc = [];
deltaDiffAcc = [];
muAcc = [];
H = H0 + Himp;

% BdG iterations
for i = 1:maxLoop
    KE = H - mu*eye(nBands);
    BdGMatrix = [KE -delta; -delta' -KE];
    [eVector eValue] = eig(BdGMatrix);
    [En,sortIndex] = sort(real(diag(eValue)));
    eVector = eVector(:,sortIndex);
    fermi = 1./(1 + exp(En/kT));
    nUpCal = (abs(eVector(1:nBands,:)).^2)*fermi;
    nDownCal = (abs(eVector((nBands + 1):end,:)).^2)*(1 - fermi);
    deltaCal = SCInteractionMatrix.*((eVector(1:nBands,:)*(((eVector((nBands + 1):end,:))').*repmat(fermi,1,nBands))));
    deltaDiff = norm(deltaCal - delta)/norm(delta);
    nDiff = abs((1/N^2)*sum(nUpCal + nDownCal) - n0)/n0;
    if (nDiff < nTol) && (deltaDiff < deltaTol)
       break % go out of loop if self-consistency is achieved
    end
    % update
    beta =  beta1 + (beta2 - beta1).*rand(1); 
    nUp = beta*nUp + (1-beta)*nUpCal;
    nDown = beta*nDown + (1-beta)*nDownCal;
    delta = beta*delta + (1-beta)*deltaCal;   
    nAvg = (1/N^2)*(sum(nUp + nDown)); 
    mu = mu - alpha*(nAvg - n0);
    nAcc = [nAcc; nAvg];
    deltaMaxAcc = [deltaMaxAcc; max(max(delta))];
    deltaMinAcc = [deltaMinAcc; min(min(delta))]; 
    muAcc = [muAcc; mu];
    deltaDiffAcc = [deltaDiffAcc; deltaDiff];
    disp([i nDiff deltaDiff]);
end
if i < maxLoop
    % save the converged result
    fileName = ['BdG_Impurity_FeSe_Toms_BS_6Dec13', '_N_', num2str(N),'_Vimp_', num2str(Vimp)];
    save(fileName,'nAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','muAcc','mu', 'deltaTol', 'nTol');
else
    disp('***********Not Converged**********')
end

% plot
figure;
subplot(2,2,1); plot(nUpAcc + nDownAcc); title('nAcc'); axis('square');
subplot(2,2,2); plot(muAcc); title('mu'); axis('square');
subplot(2,2,3); plot(deltaMaxAcc); title('deltaMax'); axis('square');
subplot(2,2,4); plot(deltaMinAcc); title('deltaMin'); axis('square');
figure; plot(deltaDiffAcc); title('Norm deltaDiff'); axis('square');



        

