function t=random_Vimp(Vimp,N,n0,nOrbitals,sublattice)
% saves impurity potential file for a dimer (2 impurities)
if nargin < 5
    sublattice=0;
end
if nargin < 4
    nOrbitals=5
end
rng('shuffle')
% select randomly n0 positions from N^2
[x,y]=meshgrid(1:N);
indices=randperm(N^2,n0);
% store the positions
imp_vec(:,1)=x(indices);
imp_vec(:,2)=y(indices);
% create the impurity matrix
imp_matr=repmat(Vimp*eye(nOrbitals, nOrbitals),1,1,n0);

%imp_vec(1,:)=[ 0 0];
filename=['random_impurity_',num2str(Vimp(1)),'_sublattice_',num2str(sublattice),'N_',num2str(N(1)),'n0_',num2str(n0),'.mat'];
% switch sublattice
%     case 1
%         % default with first iron in right upper corner (FeSe)
%         % allow for different potentials
%         if numel(Vimp)==1
%             impPotential = Vimp*eye(nOrbitals/2, nOrbitals/2);
%         else
%             impPotential = diag(Vimp);
%         end;
%         % to do: setup correct impurity potential
%         imp_matr(:,:,1) = [impPotential zeros(nOrbitals/2); impPotential zeros(nOrbitals/2)];
%     case -1
%         % not yet implemented!
%         if numel(Vimp)==1
%             impPotential = Vimp*eye(nOrbitals/2, nOrbitals/2);
%         else
%             impPotential = diag(Vimp);
%         end;
%         imp_matr(:,:,1) = [zeros(nOrbitals/2) impPotential; zeros(nOrbitals/2) impPotential];
%     case 0
%         % no sublattice (5 orbital, or 1 orbital)
%         if numel(Vimp)==1
%             imp_matr(:,:,1) = Vimp*eye(nOrbitals, nOrbitals);
%         else
%             imp_matr(:,:,1) = diag(Vimp);
%         end;        
%         imp_vec(2,:)=[ 1 0];
%         imp_matr(:,:,2) = diag(Vimp);
% end;
save(filename,'imp_vec','imp_matr');

