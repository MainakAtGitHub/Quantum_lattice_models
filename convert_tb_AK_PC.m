function [TBparameters,latticeVector]=convert_tb_AK_PC(tbfile)
tb=load(tbfile,'-ascii');
%cut=0;
nOrb=max(max(tb(:,4:5)))
sztb=size(tb);
TBparameters=[];
latticeVector=[];
for n=1:sztb(1)
    if ~isempty(TBparameters)
        n
        [ rowColIdx ] = findRowOrColumnInMatrix( latticeVector,tb(n,1:3));
    else
        rowColIdx=[];
    end;
    if isempty(rowColIdx)
        % add the lattice vector
        latticeVector=[latticeVector;tb(n,1:3)];
        rowColIdx=size(latticeVector,1);   
        TBparameters(nOrb,nOrb,rowColIdx)=0;
    end;
    TBparameters(tb(n,4),tb(n,5),rowColIdx)=TBparameters(tb(n,4),tb(n,5),rowColIdx)+tb(n,6);
end;
sublattice=1;
if nOrb==5
    sublattice=0;
elseif nOrb==1
    sublattice=0;
end;
save([tbfile,'.mat'],'TBparameters','latticeVector','sublattice');
%tb=[];
%sztb=size(TBparameters);
%for mu=1:sztb(1)
%    for nu=1:sztb(2)
%        for r=1:sztb(3)
%            if abs(TBparameters(mu,nu,r))>cut
%                tb=[tb;latticeVector(r,:),0,mu,nu,TBparameters(mu,nu,r)];
%            end;
%        end
%    end
%end