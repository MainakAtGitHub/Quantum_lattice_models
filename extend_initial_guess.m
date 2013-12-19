function f=extend_initial_guess(largeguess,smallguess,output_file,nOrbitals)

% 3. Improve inital guess
%       use impurity results for smaller system size as an initial guess
%       for larger system.
if nargin < 3
    output_file=[smallguess,'extended'];
end;
if nargin < 4
    nOrbitals=10;
end;

% input files
%load TB_hamiltonian_FeSe_2D.mat
%load Gamma_FeSe_Toms_BS_6Dec13_GammaCut_3
load(largeguess,'-mat');
%load BdG_homogeneous_FeSe_Toms_BS_6Dec13_GammaCut_3_N_13.mat
deltaH = delta; 
muH = mu;
N = sqrt(size(delta,1)/nOrbitals);
clear delta mu;
load(smallguess,'-mat');
%load BdG_Impurity_FeSe_Toms_BS_6Dec13_GammaCut_3_N_9_Vimp_0.4.mat
deltaI = delta;
muI = mu;
Nold=sqrt(size(deltaI,1)/nOrbitals);

% initial guess
mu = muI;
nBands = N^2*nOrbitals;
% maybe copy also nUp and nDown, or one has to run a convergence of the
% chemical potential only before
nUp = .6*ones(nBands,1);
nDown = .6*ones(nBands,1);
delta = deltaH;
Nchange = (N - Nold)/2;
for ix = 1:Nold
    for iy = 1:Nold
        for jx = 1:Nold
            for jy = 1:Nold
                [iRangeOld, jRangeOld] = find_lattice_translation_index(Nold, ...
                    nOrbitals, [ix iy], [jx jy]);
                [iRange, jRange] = find_lattice_translation_index(N, ...
                    nOrbitals, [ix iy] + Nchange, [jx jy] + Nchange);
                delta(iRange, jRange) = deltaI(iRangeOld, jRangeOld);
            end
        end
    end
end
% store the guess
%output_file = ['BdG_Impurity_FeSe_Toms_BS_6Dec13_GammaCut_3', '_N_', num2str(N),'_Vimp_', num2str(Vimp),'.mat'];
save(output_file,'delta','mu','nUp','nDown');
f=1;
