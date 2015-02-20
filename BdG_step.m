function [ nUpCal, nDownCal, deltaCal] = BdG_step( KE,delta, kT,nBands, SCInteractionMatrix,mKE)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if ~exist('mKE','var')
        BdGMatrix = [KE -delta; -delta' -KE];
    else
        BdGMatrix = [KE -delta; -delta' mKE];
    end;
    tic
    [eVector, eValue] = eig(BdGMatrix);
    toc
    % save some memory for following commands (here we need to save three
    % full arrays such that we get in MB:
    % 3*(2*N^2*nOrbitals)^2*8/1024/1024 (3.6G for N=25, 470M for N=15)
    %clear BdGMatrix
    [En, sortIndex] = sort(real(diag(eValue)));
    % save some memory for following commands
    %clear eValue
    eVector = eVector(:,sortIndex);
    fermi = 1./(1 + exp(En/kT));
    nUpCal = (abs(eVector(1:nBands,:)).^2)*fermi;
    nDownCal = (abs(eVector((nBands + 1):end,:)).^2)*(1 - fermi);
    global fullgamma
    if ~fullgamma
        deltaCal = SCInteractionMatrix.*((eVector(1:nBands,:)*(((eVector((nBands + 1):end,:))').*repmat(fermi,1,nBands))));
        % debuging code
        %for n=1:100
        %    dc1(n)=SCInteractionMatrix(1,n).*((eVector(n,:)*(((eVector((nBands + 1),:))').*repmat(fermi,1,1))));
        %end;
        %dc1=SCInteractionMatrix(1,:).*((eVector(1,:)*(((eVector((nBands + 1),:))').*repmat(fermi,1,nBands))));
    else
        deltaCal=zeros(nBands,nBands);
        nOrb=size(SCInteractionMatrix.int,1);
        N=int32(sqrt(nBands/nOrb));
        tic
        deltaCal=delta_full_mex(eVector,fermi,int32(nBands),int32(nOrb),int32(N),SCInteractionMatrix.int(:),int32(SCInteractionMatrix.latt),deltaCal);
        %                deltaCal=delta_full(eVector,fermi,int32(nBands),int32(nOrb),int32(N),SCInteractionMatrix.int(:),int32(SCInteractionMatrix.latt),deltaCal);
        toc
    end
end

