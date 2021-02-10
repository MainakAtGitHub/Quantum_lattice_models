function [ nUpCal, nDownCal, UpDown_deltaCal,DownUp_deltaCal,UpUp_deltaCal,DownDown_deltaCal, En, TotKE, nAnoUpDownCal, nAnoDownUpCal] = ...
    BdG_step2( KE,UpDown_delta,DownUp_delta,UpUp_delta,DownDown_delta, kT,nBands, SCInteractionMatrix,mKE,H_off_up,H_off_down,...
    SCInteractionMatrix1,SCInteractionMatrix2,SCInteractionMatrix3)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     if ~exist('mKE','var')
%         BdGMatrix = [KE, -delta; -delta', -conj(KE)];
%     else
%         BdGMatrix = [KE, H_off_up; H_off_down, mKE];
%     end;
    BdGMatrix = [[KE, H_off_up;...
        H_off_down, mKE],...
        [-UpUp_delta, -(UpDown_delta);...
        -(DownUp_delta),-DownDown_delta];...
        [-UpUp_delta, -(UpDown_delta);...
        -(DownUp_delta),-DownDown_delta]',...
        -transpose([KE, H_off_up;...
        H_off_down, mKE])];
    %tic
    %BdGMatrix1=sparse(BdGMatrix);
    %toc
    % save some memory by giving back the eigenvalues in a vector
    % does not work under matlab v8.1 or smaller
    %[eVector, eValue] = eig(BdGMatrix,'vector');
    global fullgamma cutek dress
    nargout_tmp=nargout;
    gap_equation2;
    global saveEigVec;
    if saveEigVec
        global BdGfileName
        uS = ((eVector(1:2*nBands,:)));
        vS = ((eVector((2*nBands + 1):end,:)));
        save([BdGfileName,'EigVecs'], 'uS','vS');
    end
end