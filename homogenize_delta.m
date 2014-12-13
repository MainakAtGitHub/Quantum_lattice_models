function delta_hom=homogenize_delta(delta,latticeVectorsSC,nOrb,N)
% calculate average delta over all lattice positions and put this value to
% the homogenized delta
szlV=size(latticeVectorsSC);
% initialize variable
del=zeros(nOrb,nOrb,szlV(1));
% run over all lattice translations of the pairing interaction
for n=1:szlV(1)
    deltamatrix=zeros(nOrb,nOrb);
    % add up all contributions
    for m=0:N^2-1
        ind1=mod(m+latticeVectorsSC(n,1),N^2);
        ind2=mod(m+latticeVectorsSC(n,2)*N,N^2);
        deltamatrix=deltamatrix+delta((1:nOrb)+ind1*nOrb,(1:nOrb)+ind2*nOrb);
    end;
    % take average over the lattice
    del(:,:,n)=deltamatrix/N^2;
end;
% put in the averaged gaps to the new guess
delta_hom = lattice_translation(N, del, latticeVectorsSC);
