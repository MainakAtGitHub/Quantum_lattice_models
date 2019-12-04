function f=generate_initialguess(N,Gamma_file,n0_orb,mu,random_n,neel_n)
% generate an initialguess for a system of size N times N
% (MAINAK) 
% random_n can be false or true, neel_n can be small number like
% 0.01, 0.02 etc which says how strongly neel_n initial deviation from
% half-filled nUp and nDown we want
% (MAINAK)

if nargin < 6
    neel_n=false;
end

if nargin < 5
    random_n=false;
end
% load the pairing interaction to set the superconducting order parameters
% according the expected structure
load(Gamma_file,'-mat')
% some special case for pairing interactions in 1Fe notation
if exist('delta1Fe','var')
    Gamma=delta1Fe;
    latticeVectorsSC=latticeVectors1Fe;
end
% populate a matrix according to the lattice translation (repetition of
% corresponding order parameters for a N times N system

SCInteractionMatrix = lattice_translation(N, Gamma, latticeVectorsSC);
nOrbitals=size(Gamma,1);
if nargin<3
    % no. of valence electrons per unit cell and orbital, default value for
    % FeSC, "half filling" would be at n0=1=0.5+0.5 (per spin)
    n0_orb = 1.2;
end

nBands = N^2*nOrbitals;


if exist('nup1','var')
    % if the Gamma_file also contains the fillings nup1 in the homogeneous
    % case, start putting this as initial density; yields faster convergence
    % since already converged in the case without superconductivity    nUp=repmat(nup1,N^2,1);
    nDown=repmat(ndown1,N^2,1);
    % clean values, i.e. remove values smaller than 1e-8
    SCInteractionMatrix=(abs(real(SCInteractionMatrix))>1e-8).*real(SCInteractionMatrix)+(abs(imag(SCInteractionMatrix))>1e-8).*imag(SCInteractionMatrix);
else
    % otherwise, just put constant density in each orbital
    nUp = n0_orb/2*ones(nBands,1);
    nDown =n0_orb/2*ones(nBands,1);
    
    % Mainak
    if random_n
        n_randdd = 0.01*n0_orb*rand(length(nUp),1);
        nUp = nUp+n_randdd;
        nDown = nDown-n_randdd;
    end
    % Mainak
%     % Mainak
%     if neel_n
%         sqr_nUp=zeros(N*nOrbitals);
%         sqr_nDown=zeros(N*nOrbitals);
%         for ii=1:N
%             for jj=1:N
%                 sqr_nUp(ii,(jj-1)*nOrbitals+1:jj*nOrbitals)=max(0,(-1)^(ii+jj));
%                 sqr_nDown(ii,(jj-1)*nOrbitals+1:jj*nOrbitals)=max(0,(-1)^(ii+jj+1));
%             end
%         end
%         sqr_nUp=sqr_nUp';
%         sqr_nDown=sqr_nDown';
%         nUp=n0_orb*sqr_nUp(:);
%         nDown=n0_orb*sqr_nDown(:);
%     end
%     % Mainak
    % Mainak
    if neel_n > 0
        what_frac = neel_n;
        sqr_nUp=zeros(N*nOrbitals);
        sqr_nDown=zeros(N*nOrbitals);
        for ii=1:N
            for jj=1:N
                sqr_nUp(ii,(jj-1)*nOrbitals+1:jj*nOrbitals)=(-1)^(ii+jj);
                sqr_nDown(ii,(jj-1)*nOrbitals+1:jj*nOrbitals)=(-1)^(ii+jj+1);
            end
        end
        sqr_nUp=sqr_nUp';
        sqr_nDown=sqr_nDown';
        nUp=nUp+what_frac*n0_orb*sqr_nUp(:);
        nDown=nDown+what_frac*n0_orb*sqr_nDown(:);
    end
    % Mainak

    if nargin < 4
        mu=0;
    end
    % put in a small gap of 0.05 according to the structure in the pairing
    SCInteractionMatrix=0.05*SCInteractionMatrix;
end
% use the correct variable name
delta=SCInteractionMatrix;
% clear variable to save memory
clear SCInteractionMatrix;
% save the initial guess
if size(delta,1)>11000
    % use new version for large system sizes (file limitation)
    save([Gamma_file,'guess',num2str(N)],'delta','mu','nUp','nDown','-v7.3');
    % also generate a random guess for checking that the convergence works
    % independent of initial guess
    delta=delta.*(rand(size(delta))*2-1);
    save([Gamma_file,'guess',num2str(N),'rand'],'delta','mu','nUp','nDown','-v7.3');
else
    % use default format otherwise
    save([Gamma_file,'guess',num2str(N)],'delta','mu','nUp','nDown');
    delta=delta.*(rand(size(delta))*2-1);
    save([Gamma_file,'guess',num2str(N),'rand'],'delta','mu','nUp','nDown');
end
