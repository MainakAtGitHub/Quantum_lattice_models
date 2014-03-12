function f=generate_initialguess(N,Gamma_file,n0)
load(Gamma_file,'-mat')
if exist('delta1Fe','var')
    Gamma=delta1Fe;
    latticeVectorsSC=latticeVectors1Fe;
end;
SCInteractionMatrix = lattice_translation(N, Gamma, latticeVectorsSC);
nOrbitals=size(Gamma,1);
if nargin<3
    n0 = 1.2; % no. of valence electrons per unit cell
end;
nBands = N^2*nOrbitals;
if exist('nup1','var')
    nUp=repmat(nup1,N^2,1);
    nDown=repmat(ndown1,N^2,1);
    % clean values
   SCInteractionMatrix=(abs(real(SCInteractionMatrix))>1e-8).*real(SCInteractionMatrix)+(abs(imag(SCInteractionMatrix))>1e-8).*imag(SCInteractionMatrix);
else
    nUp = n0/2*ones(nBands,1);
    nDown =n0/2*ones(nBands,1);
    mu=0;
    SCInteractionMatrix=0.05*SCInteractionMatrix;
end;
delta=SCInteractionMatrix;
save([Gamma_file,'guess',num2str(N)],'delta','mu','nUp','nDown');
delta=delta.*(rand(size(delta))*2-1);
save([Gamma_file,'guess',num2str(N),'rand'],'delta','mu','nUp','nDown');
