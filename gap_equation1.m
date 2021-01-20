  try
        tic
        [eVector, eValue] = eig(BdGMatrix);
        toc
%       clear BdGMatrix
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
    % dressing of the eVector
    if ~isempty(dress)
        eVector=eVector.*repmat(dress',nBands*2/numel(dress),nBands*2);
    end;
    eVector = eVector(:,sortIndex);
    fermi = 1./(1 + exp(En/kT));
    if ~fullgamma
       deltaCal = SCInteractionMatrix.*((eVector(1:nBands,:)*(((eVector((nBands + 1):end,:))').*repmat(fermi,1,nBands))));
        %%%%%%%%%%%%%%%%Mainak
        nAnoUpDownCal = diag(eVector(1:nBands,:)*(((eVector((nBands + 1):end,:))').*repmat(fermi,1,nBands)));
        nAnoDownUpCal = conj(nAnoUpDownCal);
%         nAnoDownUp = (eVector(1:nBands,:))'.*eVector((nBands + 1):end,:).*repmat(fermi,1,nBands);
%         if sum(abs(nAnoUpDown-conj(nAnoDownUp))) > 1e-5
%             disp('Nonhermitian hamiltonian!');
%         end
        %%%%%%%%%%%%%%%%%%%%Mainak
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
        %fermi=-.5+fermi;
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
    nUpCal = (abs(eVector(1:nBands,:)).^2)*fermi;
    nDownCal = (abs(eVector((nBands + 1):end,:)).^2)*(fermi);
    %anomalousncalc
    % problem in matlab R2018b and later
    if nargout_tmp > 6 
    TotKE = sum(diag((eVector)'*BdGMatrix*(eVector)).*fermi);%.*tanh(En/kT));
    end
    % Mainak
%     uCal = (abs(eVector(1:nBands,floor(nBands/2)).^2));%*ones(size(fermi));
%     vCal = (abs(eVector((nBands + 1):end,floor(nBands/2))).^2);%*ones(size(1-fermi));
    % Mainak
