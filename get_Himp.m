function Himp=get_Himp(Vimp,N,nOrbitals,sublattice,randompot,BdGfileName,dislocation_length)
if nargin<7
    dislocation_length=false;
end
if nargin <5
    randompot=0;
end;

% get_Himp(Vimp,N,nOrbitals,sublattice) gives back the Hamiltonian for the
% impurity or sets of impurities given in the input file Vimp

% Mainak%%% temporarily changed due to disloc
if numel(N)==1
    N=[N,N];
end
if true %dislocation_length>0
    nBands = (N(1)*N(2)-dislocation_length)*nOrbitals;
else
    % Mainak
    
    nBands = N(1)*N(2)*nOrbitals;
end
Himp = zeros(nBands);
% rectangular by default
if numel(N)==1
     N=[N,N];
end
impCell = [ceil(N(1)/2) ceil(N(2)/2)];
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
        % does the same as the following lines
        % cellvector=periodic_latticevectors(cellvector,N);
        cellvector(1)=mod(cellvector(1)-1,N(1))+1;
        cellvector(2)=mod(cellvector(2)-1,N(2))+1;
        [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, cellvector, cellvector);
        if ~iscell(imp_matr)
            % only on-site potentials
            if size(imp_matr,1)<iRange
                % put a 5 orbital potential on the first lattice position
                iRange=iRange(1:size(imp_matr,1));
                jRange=jRange(1:size(imp_matr,2));
                Himp(iRange, jRange)=Himp(iRange, jRange)+imp_matr(:,:,n);
            else
                Himp(iRange, jRange)=Himp(iRange, jRange)+imp_matr(:,:,n);
            end;
        else
            % read impurity potential given in filename (slow, because
            % multiple times reading the same input, but working for now)
            load(imp_matr{n});
            % add the corresponding impurity hoppings, potentials
            Himp=Himp+lattice_single(N, Impparameters, ImpVector,cellvector);
        end;
    end;
end;
if randompot>0
    % add a random potential for splitting of degenerate states
    % try to load previously used potential information
    try
        [pathstr,name,ext] =fileparts(BdGfileName);
        potentialfile=[pathstr,filesep,name,'rnd',ext];
        load(potentialfile);
        Hrand=0*Himp;
        for n=1:size(randomcells,1)
            [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, randomcells(n,:),randomcells(n,:));
            Hrand(iRange, jRange)=diag(randompotentials(n,:));
        end;
    catch
        fraction=2;
        threshold=0.5;
        disp('Generating random potentials on every ',num2str(fraction),' site with threshold of ',num2str(threshold));
        Hrand=0*Himp;
        randomcells=[];
        randompotentials=[];
        for nx=1:N(1)
            for ny=1:N(2)
                % decide whether to put a random potential
                rnd=rand(1);
                if rnd >1/fraction
                    randomcells=[randomcells;nx,ny];
                    [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, [nx,ny], [nx,ny]);
                    rndpot=rand(1,nOrbitals)-0.5;
                    rndpot=rndpot.*(abs(rndpot)>0.5*threshold);
                    randompotentials=[randompotentials;rndpot];
                    Hrand(iRange, jRange)=diag(rndpot);
                end;
            end
        end;
        save(potentialfile,'randompotentials','randomcells');
    end;
    Himp=Himp+Hrand*randompot;
end;




