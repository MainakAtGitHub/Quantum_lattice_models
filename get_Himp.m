function Himp=get_Himp(Vimp,N,nOrbitals,sublattice)
% get_Himp(Vimp,N,nOrbitals,sublattice) gives back the Hamiltonian for the
% impurity or sets of impurities given in the input file Vimp
nBands = N^2*nOrbitals;
Himp = zeros(nBands);
impCell = [ceil(N/2) ceil(N/2)];
% allow for general impurity potentials
if ~ischar(Vimp)
    [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, impCell, impCell);
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
            Himp(iRange, jRange) = [impPotential zeros(nOrbitals/2); zeros(nOrbitals/2) zeros(nOrbitals/2)];
        case -1
            if numel(Vimp)==1
                impPotential = Vimp*eye(nOrbitals/2, nOrbitals/2);
            else
                impPotential = diag(Vimp);
            end;
            Himp(iRange, jRange) = [zeros(nOrbitals/2) zeros(nOrbitals/2); zeros(nOrbitals/2) impPotential];
            % to be implemented
        case 0
            % no sublattice (5 orbital)
            if numel(Vimp)==1
                Himp(iRange, jRange) = Vimp*eye(nOrbitals, nOrbitals);
            else
                Himp(iRange, jRange) = diag(Vimp);
            end;
    end;
else
    load(Vimp,'-mat')
    % now we have a set of impurity matrices imp_matr
    % together with some lattice vectors imp_vec
    % here we set up Himp directly
    numimp=size(imp_vec,1);
    for n=1:numimp
        cellvector=impCell+imp_vec(n,:);
        cellvector(1)=mod(cellvector(1)-1,N)+1;
        cellvector(2)=mod(cellvector(2)-1,N)+1;
        [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, cellvector, cellvector);
        Himp(iRange, jRange)=imp_matr(:,:,n);
    end;
end;



