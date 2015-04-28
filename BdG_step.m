function [ nUpCal, nDownCal, deltaCal] = BdG_step( KE,delta, kT,nBands, SCInteractionMatrix,mKE)
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
    try
        tic
        [eVector, eValue] = eig(BdGMatrix);
        toc
        clear BdGMatrix
        [En, sortIndex] = sort(diag(eValue));
    catch
        disp('error with buildin eig, using eigen3 library now')
        clear eVector eValue
        tic
        [eVector, eValue] =    SelfAdjointEigenSolver(BdGMatrix);
        toc
        [En, sortIndex] = sort(eValue);
    end

       % sparse matrix eigenvalue calculation is significantly slower!
        %[eVector, eValue] = eigs(BdGMatrix1,nBands);
    % save some memory for following commands (here we need to save three
    % full arrays such that we get in MB:
    % 3*(2*N^2*nOrbitals)^2*8/1024/1024 (3.6G for N=25, 470M for N=15)
    % save some memory for following commands
    clear eValue
    eVector = eVector(:,sortIndex);
    fermi = 1./(1 + exp(En/kT));
    nUpCal = (abs(eVector(1:nBands,:)).^2)*fermi;
    nDownCal = (abs(eVector((nBands + 1):end,:)).^2)*(1 - fermi);
    global fullgamma cutek
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
        % symmetrize with respect to particles/antiparticles
        fermi=-.5+fermi;
        tic
        if ~isnan(cutek)
            % restrict summation over finite range of energies
            maxdelta=max(max(abs(delta)));
            if cutek < 0
                cut=-maxdelta*cutek;
            else
                cut=cutek;
                if cut < 24*maxdelta
                    disp(['largest gap is ',num2str(maxdelta),' cut for summation set too small, setting to',num2str(24*maxdelta)]);
                    cut=24*maxdelta;
                end;
            end;
            eklist=find(abs(En)<cut);
            eVector=eVector(:,eklist);
            fermi=fermi(eklist);
            disp(['summing over ',num2str(numel(fermi)), ' instead of ',num2str(2*nBands),' energies']);
        end;
        deltaCal=delta_full_mex(eVector,fermi,int32(nBands),int32(nOrb),int32(N),SCInteractionMatrix.int(:)+0*1i,int32(SCInteractionMatrix.latt),deltaCal);
        %deltaCal=delta_full(eVector,fermi,int32(nBands),int32(nOrb),int32(N),SCInteractionMatrix.int(:),int32(SCInteractionMatrix.latt),deltaCal);
        toc
    end
end

