function [ nUpCal_super, nDownCal_super, deltaCal_super, En] = BdG_step_super( HSuper,delta, kT,nBands, SCInteractionMatrix,BZ,HImpurity,mu,mKE)
% supercell version of BdG_step
% BZ contains information about the number of k-points to be used and
% symmetries to be explored
global fullgamma cutek
%   Detailed explanation goes here
startindex=1;
endindex=size(BZ.k,1);
% repeated code here
if fullgamma
	maxHop = max(max(abs(SCInteractionMatrix.latt)));
    nOrbitals=size(SCInteractionMatrix.int,1);
else
    disp('not implemented');
    %maxHop = max(max(abs(Gammafull.latt)));
end
[deltaSuper,superLatticeVectors] = supercell_delta(nOrbitals, delta, maxHop);
deltaCal_super=zeros(nBands,nBands);
nUpCal_super=0;
nDownCal_super=0;
for index=startindex:endindex
    if ~exist('mKE','var')
       % iKy= mod(index-1,BZ.M)+1;
       % iKx= ceil(index/BZ.M);
        nSuperCells = size(superLatticeVectors,1);
        k = BZ.k(index,:);
        disp(k);
        kSpaceHopping = 0;
        kSpaceHoppingc=0;
        kSpaceGap = 0;
        for iUnitCell = 1:nSuperCells
            iLatticeVector = superLatticeVectors(iUnitCell,:);
            kSpaceHopping = kSpaceHopping + HSuper(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
            kSpaceHoppingc = kSpaceHoppingc + HSuper(:,:,iUnitCell)*exp(-1i*(iLatticeVector*k'));
            kSpaceGap = kSpaceGap + deltaSuper(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
        end
        KESuper = kSpaceHopping + HImpurity-mu*eye(nBands);
        KESuperc = kSpaceHoppingc + HImpurity;
        KESuper=0.5*[KESuper+KESuper'];
        KESuperc=0.5*[KESuperc+KESuperc'];
       % BdGMatrix = [KE, -delta; -delta', -conj(KE)];
        BdGMatrix = [KESuper -kSpaceGap; -kSpaceGap' -conj(KESuperc)];
    else
        disp('not implemented');
        %BdGMatrix = [KE, -delta; -delta', mKE];
    end;
    %tic
    %BdGMatrix1=sparse(BdGMatrix);
    %toc
    % save some memory by giving back the eigenvalues in a vector
    % does not work under matlab v8.1 or smaller
    %[eVector, eValue] = eig(BdGMatrix,'vector');
    gap_equation;
    deltaCal_super=deltaCal_super+BZ.weight(index)*deltaCal;
    nUpCal_super=nUpCal_super+BZ.weight(index)*nUpCal;
    nDownCal_super=nDownCal_super+BZ.weight(index)*nDownCal;
end