function f=generate_initialguess(N,Gamma_file)
load(Gamma_file,'-mat')
SCInteractionMatrix = lattice_translation(N, Gamma, latticeVectorsSC);
nOrbitals=size(Gamma,1);
    n0 = 1.2; % no. of valence electrons per unit cell
nBands = N^2*nOrbitals;
nUp = n0/2*ones(nBands,1);
nDown =n0/2*ones(nBands,1);
mu=0;
delta=0.05*SCInteractionMatrix;
save([Gamma_file,'guess',num2str(N)],'delta','mu','nUp','nDown');
