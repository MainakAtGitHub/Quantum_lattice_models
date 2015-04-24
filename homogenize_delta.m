function delta_hom=homogenize_delta(delta,latticeVectorsSC,nOrb,N)
% calculate average delta over all lattice positions and put this value to
% the homogenized delta
szlV=size(latticeVectorsSC);
% initialize variable
del=zeros(nOrb,nOrb,szlV(1));
singlet=false;
% run over all lattice translations of the pairing interaction
for n=1:szlV(1)
    deltamatrix=zeros(nOrb,nOrb);
    % add up all contributions
    for m=0:N^2-1
        x=fix(m/N);
        y=mod(m,N);
        x1=x+latticeVectorsSC(n,1);
        y1=y+latticeVectorsSC(n,2);
        x1=mod(x1+floor(y/N),N);
        y1=mod(y1,N);
        ind2=mod(m,N^2)*nOrb;
        ind1=mod(x1*N+y1,N^2)*nOrb;
        %tmp1=mod(fix(m+latticeVectorsSC(n,2)/N),N)*0;
        %ind1=mod(m+(latticeVectorsSC(n,1)+tmp1)*N+latticeVectorsSC(n,2),N^2)*nOrb;
        deltamatrix=deltamatrix+delta((1:nOrb)+ind1,(1:nOrb)+ind2);
        %disp([ind1,ind2,delta(ind1+1,ind2+1)]);
        if singlet
            % project to singlet channel (not so easy, we need to consider
            % also the properties of the orbitals under r -> -r
            ind1m=mod(m-latticeVectorsSC(n,1)*N,N^2);
            ind2m=mod(m-latticeVectorsSC(n,2),N^2);
            deltamatrix=deltamatrix+delta((1:nOrb)+ind1m*nOrb,(1:nOrb)+ind2m*nOrb);
        end;
    end;
    % take average over the lattice
    del(:,:,n)=deltamatrix/N^2;
end;
% put in the averaged gaps to the new guess
delta_hom = lattice_translation(N, del, latticeVectorsSC);
if singlet
    delta_hom=0.5*delta_hom;
end;