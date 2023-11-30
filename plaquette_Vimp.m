function t=plaquette_Vimp(Vimp,nOrbitals,sublattice)
% saves impurity potential file for a plaquette (4 impurities)
if nargin < 3
    sublattice=1
end;
if nargin < 2
    nOrbitals=10
end;
imp_vec(1,:)=[ 0 0];
filename=['plaquette_impurity_',num2str(Vimp(1)),'_sublattice_',num2str(sublattice),'.mat']
switch sublattice
    case 1
        % default with first iron in right upper corner (FeSe)
        % allow for different potentials
        if numel(Vimp)==1
            impPotential = Vimp*eye(nOrbitals/2, nOrbitals/2);
        else
            impPotential = diag(Vimp);
        end;
        % to do: setup correct impurity potential
        imp_matr(:,:,1) = [impPotential zeros(nOrbitals/2); zeros(nOrbitals/2) impPotential];
        imp_vec(2,:)=[0 1];
        imp_matr(:,:,2) = [impPotential zeros(nOrbitals/2); zeros(nOrbitals/2) zeros(nOrbitals/2)];
        imp_vec(3,:)=[1 0];
        imp_matr(:,:,3) = [0*impPotential zeros(nOrbitals/2);  zeros(nOrbitals/2) impPotential];
    case -1
        % not yet implemented!
        disp('Error, not yet implemented');
        if numel(Vimp)==1
            impPotential = Vimp*eye(nOrbitals/2, nOrbitals/2);
        else
            impPotential = diag(Vimp);
        end;
        imp_matr = [impPotential zeros(nOrbitals/2); zeros(nOrbitals/2) impPotential];
        % to be implemented
    case 0
        % no sublattice (5 orbital, or 1 orbital)
        if numel(Vimp)==1
            imp_matr(:,:,1) = Vimp*eye(nOrbitals, nOrbitals);
        else
            imp_matr(:,:,1) = diag(Vimp);
        end;        
        imp_vec(2,:)=[ 1 0];
        imp_matr(:,:,2) = diag(Vimp);
        imp_vec(3,:)=[ 0 1];
        imp_matr(:,:,3) = diag(Vimp);  
        imp_vec(4,:)=[ 1 1];
        imp_matr(:,:,4) = diag(Vimp);
end;
save(filename,'imp_vec','imp_matr');

