function t=trivial_Vimp(Vimp,nOrbitals,sublattice)
% get_Himp(Vimp,N,nOrbitals,sublattice) gives back the Hamiltonian for the
% impurity or sets of impurities given in the input file Vimp
if nargin < 3
    sublattice=1
end;
if nargin < 2
    nOrbitals=10
end;
imp_vec=[ 0 0];
filename=['trivial_impurity_',num2str(Vimp(1)),'_sublattice_',num2str(sublattice),'.mat']
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
        imp_matr = [impPotential zeros(nOrbitals/2); zeros(nOrbitals/2) zeros(nOrbitals/2)];
    case -1
        if numel(Vimp)==1
            impPotential = Vimp*eye(nOrbitals/2, nOrbitals/2);
        else
            impPotential = diag(Vimp);
        end;
        imp_matr = [zeros(nOrbitals/2) zeros(nOrbitals/2); zeros(nOrbitals/2) impPotential];
        % to be implemented
    case 0
        % no sublattice (5 orbital)
        if numel(Vimp)==1
            imp_matr = Vimp*eye(nOrbitals, nOrbitals);
        else
            imp_matr = diag(Vimp);
        end;
end;
save(filename,'imp_vec','imp_matr');

