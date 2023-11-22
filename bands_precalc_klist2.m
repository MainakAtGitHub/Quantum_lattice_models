function [kSpaceEigenValuesNormal,kSpaceEigenVectorsNormal]=bands_precalc_klist2(inputfile,klist,onsite_SOC,switch_on_nm_dg_br)

if nargin < 4
    switch_on_nm_dg_br = 0;
end
if nargin < 3 
    onsite_SOC=0; % onsite_SOC used as boolean
end
% for given inputfile (*.txt) and an array containing the momentum space Hamiltonian, calculate the eigenvalues and eigenvectors
% of the band structure

% call the function that does the Fourier transform
[kSpaceHamiltonian]=Hk_klist2(inputfile,klist);

% % % % % % % by default set the chemical potential to zero
% % % % % %  mu=0;%mu=1.062;
% % % % % %     % load the input file to set the variables, gave up the old .mat file
% % % % % %     % format
     read_input_file=inputfile;
     read_input;
     read_input_file
% % % % % % 
% % % % % % 
% % % % % % % load the hoppings from the file
% % % % % % load(TB_file,'-mat');
% % % % % % nOrbitals = size(TBparameters,1);
% % % % % % 
% % % % % % % number of hoppings in band structure
% % % % % % nUnitCells = size(latticeVector,1);
% % % % % % 
% % % % % % % allow to set the chemical potential
% % % % % % load(BdGfileName,'-mat');  %%%%%remember to read BdGfileName for mu
% % % % % % mu=real(mu);
% % % % % % TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) = ...
% % % % % % TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) - mu*eye(nOrbitals);
% % % % % % 
% % % % % % % fix for the lattice vectors (no 3rd component allowed)
% % % % % % latticeVector=latticeVector(:,1:2);

% check the size of the input array
szk=size(kSpaceHamiltonian);
if ~(szk(2)==nOrbitals)
    disp('mismatch of number of orbitals, check input')
end
% initialize the result variables
kSpaceEigenValuesNormal = zeros(szk(1), nOrbitals*(onsite_SOC+1));
kSpaceEigenVectorsNormal = zeros(szk(1), nOrbitals*(onsite_SOC+1), nOrbitals*(onsite_SOC+1));

%iterate over the momentum vectors
if ~onsite_SOC
    for iKx = 1:szk(1)
        % diagonalize the Hamiltonian
        [eigVectorNormal, eigValueNormal] = eig(squeeze(kSpaceHamiltonian(iKx,:,:)));
        % sort the eigenvalues and give back the order
        [eigValueKNormal, sortingIndexNormal] = sort(real(diag(eigValueNormal)));
        % put eigenvectors in the corresponding order
        eigVectorKNormal = (eigVectorNormal(:,sortingIndexNormal))';
        % put them into the restult array
        kSpaceEigenValuesNormal(iKx, :) = eigValueKNormal;
        kSpaceEigenVectorsNormal(iKx, :,:) =  eigVectorKNormal;
    end
end
if onsite_SOC
    lmbda_SOC_onsite=spin_and_nambu_SOC;
    
    if ~exist('Orb_seq','var')
        % z^2 (1), xz (2), yz (3), x^2-y^2 (4), xy (5)
        % Orb_seq=[2,3,1,5,4];
        %Orb_seq=[2,3,4,5,1];
        %Orb_seq=[3,2,1,5,4];
        %Orb_seq=[3,2,4,5,1];
        disp('Define orbital sequence for SO coupling')
    end
    
    
    Lx5 = [[0             0  sqrt(3)*1i      0    0];...
        [0             0           0      0   1i];...
        [-sqrt(3)*1i   0           0    -1i    0];...
        [0             0          1i      0    0];...
        [0           -1i           0      0    0]];
    
    Ly5 = [[0    sqrt(3)*1i           0      0    0];...
        [-sqrt(3)*1i   0           0     1i    0];...
        [0             0           0      0   1i];...
        [0           -1i           0      0    0];...
        [0             0         -1i      0    0]];
    
    Lz5 = [[0             0           0      0    0];...
        [0             0         -1i      0    0];...
        [0            1i           0      0    0];...
        [0             0           0      0  -2i];...
        [0             0           0     2i    0]];
    
    Lx5=Lx5(:,Orb_seq); Lx5= Lx5(Orb_seq,:);
    Ly5=Ly5(:,Orb_seq); Ly5= Ly5(Orb_seq,:);
    Lz5=Lz5(:,Orb_seq); Lz5= Lz5(Orb_seq,:);
    
    HSOC_UU = Lz5(1:nOrbitals,1:nOrbitals)/2;
    HSOC_DD = -Lz5(1:nOrbitals,1:nOrbitals)/2;
    HSOC_UD = (Lx5(1:nOrbitals,1:nOrbitals)+1i*Ly5(1:nOrbitals,1:nOrbitals))/2;
    HSOC_DU = (Lx5(1:nOrbitals,1:nOrbitals)-1i*Ly5(1:nOrbitals,1:nOrbitals))/2;

    if switch_on_nm_dg_br
        nm_dg_br=10^-5; % numerical_degeneracy_breaker
    else 
        nm_dg_br=0;
    end
    disp('numerical_degeneracy_breaker implemented');
    for iKx = 1:szk(1)
        % diagonalize the Hamiltonian
        [eigVectorNormal, eigValueNormal] = eig([squeeze(kSpaceHamiltonian(iKx,:,:))+lmbda_SOC_onsite*HSOC_UU-nm_dg_br*eye(size(HSOC_UU)),                 lmbda_SOC_onsite*HSOC_UD;...
                                                              lmbda_SOC_onsite*HSOC_DU,                        squeeze(conj(kSpaceHamiltonian(iKx,:,:)))+lmbda_SOC_onsite*HSOC_DD+nm_dg_br*eye(size(HSOC_DD))]);
        % sort the eigenvalues and give back the order
        [eigValueKNormal, sortingIndexNormal] = sort(real(diag(eigValueNormal)));
        % put eigenvectors in the corresponding order
        eigVectorKNormal = (eigVectorNormal(:,sortingIndexNormal))';
        % put them into the restult array
        kSpaceEigenValuesNormal(iKx, :) = eigValueKNormal;
        kSpaceEigenVectorsNormal(iKx, :,:) =  eigVectorKNormal;
    end    
end


% % % % % % kSpaceEigenValuesNormal=kSpaceEigenValuesNormal-mu;