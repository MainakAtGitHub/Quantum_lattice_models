function [ nUpCal, nDownCal, deltaCal, En, TotKE,nAnoUpDownCal,nAnoDownUpCal] = BdG_step1( KE,delta, kT,nBands, SCInteractionMatrix,mKE,H_off_up,H_off_down)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if ~exist('mKE','var')
        BdGMatrix = [KE, -delta; -delta', -conj(KE)];
    else
        BdGMatrix = [KE, H_off_up; H_off_down, mKE];
    end;
    %tic
    %BdGMatrix1=sparse(BdGMatrix);
    %toc
    % save some memory by giving back the eigenvalues in a vector
    % does not work under matlab v8.1 or smaller
    %[eVector, eValue] = eig(BdGMatrix,'vector');
    global fullgamma cutek dress
    nargout_tmp=nargout;
    gap_equation1;
    global saveEigVec;
    if saveEigVec
        global BdGfileName
        uS = ((eVector(1:nBands,:)));
        vS = ((eVector((nBands + 1):end,:)));
        save([BdGfileName,'EigVecs'], 'uS','vS');
    end
end