function t=trivial_gamma(nOrb,neighbor,magnitude)
% generation of a trivial interaction for nOrb Orbitals
% 10 orbitals (as in FeSe)
if nargin < 2
    % for the NN bonds
    neighbor=1
end;
if nargin < 3
    magnitude=0.3
end;
if nargin < 1
    nOrb=10;
end;
switch neighbor
    case 1
        latticeVectorsSC(1,:)=[ 0 0];
        oneblock=diag(ones(1,nOrb/2));
        zblock=0*oneblock;
        Gamma(:,:,1)=[oneblock,zblock;zblock,oneblock];
        Gamma=Gamma*0.01;
        latticeVectorsSC(1,:)=[ 0 0];
        latticeVectorsSC(2,:)=[ 1 0];
        latticeVectorsSC(3,:)=[ 0 1];
        latticeVectorsSC(4,:)=[ -1 0];
        latticeVectorsSC(5,:)=[ 0 -1];
        latticeVectorsSC(6,:)=[ -1 1];
        latticeVectorsSC(7,:)=[ 1 -1];
        oneblock=diag(ones(1,5));
        zblock=0*oneblock;
        TBparameters(:,:,1)=[oneblock,oneblock;oneblock,oneblock];
        TBparameters(:,:,2)=[zblock,oneblock;zblock,zblock];
        TBparameters(:,:,3)=[zblock,zblock;oneblock,zblock];
        TBparameters(:,:,4)=[zblock,zblock;oneblock,zblock];
        TBparameters(:,:,5)=[zblock,oneblock;zblock,zblock];
        TBparameters(:,:,6)=[zblock,zblock;oneblock,zblock];
        TBparameters(:,:,7)=[zblock,oneblock;zblock,zblock];
        Gamma=TBparameters*.5;
        outputfile=['Gamma_trivial',num2str(nOrb),'.mat']
    case 2
        % NNN bond for 2 sublattice description as in FeSe (possibly also
        % correct with the switched positions)
        latticeVectorsSC(1,:)=[ 1 0];
        latticeVectorsSC(2,:)=[ 0 1];
        latticeVectorsSC(3,:)=[ -1 0];
        latticeVectorsSC(4,:)=[ 0 -1];
        Gamma(:,:,1)=eye(nOrb);
        Gamma(:,:,2)=eye(nOrb);
        Gamma(:,:,3)=eye(nOrb);
        Gamma(:,:,4)=eye(nOrb);
        % set gamma to -gamma (since there is no additional sign in the BdG
        % code)
        Gamma=magnitude*Gamma;        
        outputfile=['Gamma_trivial_NNN',num2str(nOrb),'_',num2str(magnitude),'.mat']
end
save(outputfile,'Gamma','latticeVectorsSC','-mat');