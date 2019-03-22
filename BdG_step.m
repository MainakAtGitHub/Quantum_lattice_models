function [ nUpCal, nDownCal, deltaCal, En] = BdG_step( KE,delta, kT,nBands, SCInteractionMatrix,mKE)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if ~exist('mKE','var')
        BdGMatrix = [KE, -delta; -delta', -conj(KE)];
    else
        BdGMatrix = [KE, -delta; -delta', mKE];
    end;
    %tic
    %BdGMatrix1=sparse(BdGMatrix);
    %toc
    % save some memory by giving back the eigenvalues in a vector
    % does not work under matlab v8.1 or smaller
    %[eVector, eValue] = eig(BdGMatrix,'vector');
    global fullgamma cutek dress
    gap_equation;
end

