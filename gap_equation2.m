    try
        tic
        [eVector, eValue] = eig(BdGMatrix);
        toc
%       clear BdGMatrix
%         [En, sortIndex] = sort(diag(eValue));
        En=diag(eValue);
    catch
        disp('error with buildin eig, using eigen3 library now')
        clear eVector eValue
        tic
        [eVector, eValue] = SelfAdjointEigenSolver(BdGMatrix);
        toc
%         [En, sortIndex] = sort(eValue);
        En=diag(eValue);
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
        eVector=eVector.*repmat(dress', nBands*2/numel(dress),nBands*2);
    end;
%     eVector = eVector(:,sortIndex);
    fermi = 1./(1 + exp(En/kT));
    if ~fullgamma
        UpDownExpec=((eVector(1:nBands,:)*(((eVector((3*nBands + 1):end,:))').*repmat(1-fermi,1,nBands))));
         
%%%%%%% Lines below are also crrect ways of calculating the expectation values        
%         DownUpExpec1=(conj(eVector((3*nBands + 1):end,:)*((transpose(eVector(1:nBands,:))).*repmat(fermi,1,nBands))));
%         UpDownExpec1=(conj(eVector((2*nBands + 1):3*nBands,:))*((transpose(eVector(1*nBands+1:2*nBands,:))).*repmat(fermi,1,nBands)));
        
        DownUpExpec=((eVector(1*nBands+1:2*nBands,:)*(((eVector((2*nBands + 1):3*nBands,:))').*repmat(1-fermi,1,nBands))));
        UpUpExpec=((eVector(1:nBands,:)*(((eVector((2*nBands + 1):3*nBands,:))').*repmat((1-fermi),1,nBands))));
        DownDownExpec=((eVector(nBands+1:2*nBands,:)*(((eVector((3*nBands + 1):end,:))').*repmat((1-fermi),1,nBands))));
        
        
        %%%%%%% Lines below are also crrect ways of calculating the expectation values  
%         UpUpExpec1=(((conj(eVector((2*nBands + 1):3*nBands,:)))*((transpose(eVector(1:nBands,:))).*repmat((fermi),1,nBands))));
%         DownDownExpec1=(((conj(eVector((3*nBands + 1):4*nBands,:)))*((transpose(eVector(1*nBands+1:2*nBands,:))).*repmat((fermi),1,nBands))));        
        
        ExpMatNorUpUpCal=((eVector(1:nBands,:)*(((eVector(1:nBands,:))').*repmat(fermi,1,nBands))));
        ExpMatNorUpDownCal=eVector(1*nBands+1:2*nBands,:)*(((eVector((0*nBands + 1):1*nBands,:))').*repmat(fermi,1,nBands)); % opposite convention of nAnoDownUpCal followed here
        ExpMatNorDownUpCal=eVector(0*nBands+1:1*nBands,:)*(((eVector((nBands + 1):2*nBands,:))').*repmat(fermi,1,nBands)); % opposite convention of nAnoUpDownCal followed here
        ExpMatNorDownDownCal=((eVector(nBands+1:2*nBands,:)*(((eVector(nBands+1:2*nBands,:))').*repmat(fermi,1,nBands))));
        
        %%%%test variables
%         testUpDownExpec=((eVector(1:nBands,:)*(((eVector((3*nBands + 1):end,:))').*repmat(fermi,1,nBands))));
%         testDownUpExpec=((eVector(1*nBands+1:2*nBands,:)*(((eVector((2*nBands + 1):3*nBands,:))').*repmat(fermi,1,nBands))));
%         testUpUpExpec=((eVector(1:nBands,:)*(((eVector((2*nBands + 1):3*nBands,:))').*repmat((1-fermi),1,nBands))));
%         testDownDownExpec=((eVector(nBands+1:2*nBands,:)*(((eVector((3*nBands + 1):end,:))').*repmat((1-fermi),1,nBands))));
        %%%%test variables
        
        UpDown_deltaCal=-(SCInteractionMatrix.*(UpDownExpec - DownUpExpec)*0.5);
        DownUp_deltaCal=transpose(SCInteractionMatrix.*(UpDownExpec - DownUpExpec)*0.5);
%         UpUp_deltaCal=zeros(size(UpDown_deltaCal));
%         DownDown_deltaCal=zeros(size(UpDown_deltaCal));

        
        UpUp_deltaCal1=transpose(SCInteractionMatrix1.*(UpUpExpec - DownDownExpec)*0.5);%zeros(size(UpDown_deltaCal1));
        DownDown_deltaCal1=-transpose(SCInteractionMatrix1.*(UpUpExpec - DownDownExpec)*0.5);%zeros(size(UpDown_deltaCal1));
%         UpDown_deltaCal1=zeros(size(UpUp_deltaCal1));
%         DownUp_deltaCal1=zeros(size(UpUp_deltaCal1));

        UpUp_deltaCal2=transpose(SCInteractionMatrix2.*(UpUpExpec + DownDownExpec)*0.5);%zeros(size(UpDown_deltaCal1));
        DownDown_deltaCal2=transpose(SCInteractionMatrix2.*(UpUpExpec + DownDownExpec)*0.5);%zeros(size(UpDown_deltaCal1));
%         UpDown_deltaCal2=zeros(size(UpUp_deltaCal2));
%         DownUp_deltaCal2=zeros(size(UpUp_deltaCal2));
        
        UpDown_deltaCal3=-(SCInteractionMatrix3.*(UpDownExpec + DownUpExpec)*0.5);
        DownUp_deltaCal3=transpose(SCInteractionMatrix3.*(UpDownExpec + DownUpExpec)*0.5);
%         UpUp_deltaCal3=zeros(size(UpDown_deltaCal3));
%         DownDown_deltaCal3=zeros(size(UpDown_deltaCal3));
        
        
        UpDown_deltaCal = UpDown_deltaCal + UpDown_deltaCal3;
        DownUp_deltaCal = DownUp_deltaCal + DownUp_deltaCal3;
        UpUp_deltaCal =   UpUp_deltaCal1 + UpUp_deltaCal2;
        DownDown_deltaCal = DownDown_deltaCal1 + DownDown_deltaCal2;
        
        global saveEnTot_spin_and_nambu
        if saveEnTot_spin_and_nambu
            global E_Sup
            E_Sup = sum(sum((SCInteractionMatrix.*(abs(UpDownExpec - DownUpExpec)).^2*0.5))) + sum(sum(transpose(SCInteractionMatrix.*(abs(UpDownExpec - DownUpExpec)).^2*0.5)))...
                + sum(sum(transpose(SCInteractionMatrix1.*(abs(UpUpExpec - DownDownExpec)).^2*0.5))) + sum(sum(transpose(SCInteractionMatrix1.*(abs(UpUpExpec - DownDownExpec)).^2*0.5)))...
                + sum(sum(transpose(SCInteractionMatrix2.*(abs(UpUpExpec + DownDownExpec)).^2*0.5))) + sum(sum(transpose(SCInteractionMatrix2.*(abs(UpUpExpec + DownDownExpec)).^2*0.5)))...
                + sum(sum((SCInteractionMatrix3.*(abs(UpDownExpec + DownUpExpec)).^2*0.5))) + sum(sum(transpose(SCInteractionMatrix3.*(abs(UpDownExpec + DownUpExpec)).^2*0.5)));
        end
        
%         SCInteractionMatrix.*((eVector(1:nBands,:)*(((eVector((nBands + 1):end,:))').*repmat(fermi,1,nBands))));
%         UpDown_deltaCal = SCInteractionMatrix.*((eVector(1:nBands,:)*(((eVector((3*nBands + 1):end,:))').*repmat(fermi,1,nBands))));
%        DownUp_deltaCal = SCInteractionMatrix.*((eVector(1*nBands+1:2*nBands,:)*(((eVector((2*nBands + 1):3*nBands,:))').*repmat(fermi,1,nBands))));
%        UpUp_deltaCal = 0*SCInteractionMatrix.*((eVector(1:nBands,:)*(((eVector((2*nBands + 1):3*nBands,:))').*repmat((1-fermi),1,nBands))));
%        DownDown_deltaCal = 0*SCInteractionMatrix.*((eVector(nBands+1:2*nBands,:)*(((eVector((3*nBands + 1):end,:))').*repmat((1-fermi),1,nBands))));
        %%%%%%%%%%%%%%%%Mainak
        nAnoUpDownCal=diag(ExpMatNorDownUpCal); %%%%%%TESTING aug2021
        nAnoDownUpCal=diag(ExpMatNorUpDownCal); %%%%%%TESTING aug2021
% % % %         nAnoUpDownCal = diag(eVector(0*nBands+1:1*nBands,:)*(((eVector((nBands + 1):2*nBands,:))').*repmat(fermi,1,nBands)));
% % % %         nAnoDownUpCal = diag(eVector(1*nBands+1:2*nBands,:)*(((eVector((0*nBands + 1):1*nBands,:))').*repmat(fermi,1,nBands)));%conj(nAnoUpDownCal);%diag(eVector(3*nBands+1:4*nBands,:)*(((eVector((0*nBands + 1):1*nBands,:))').*repmat(fermi,1,nBands)));
%         nAnoDownUpCal = conj(nAnoUpDownCal);
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
    nUpCal = diag(ExpMatNorUpUpCal);            %%%%%%TESTING aug2021
    if max(abs(imag(nUpCal))) > 10^-15        
        disp(['imaginary densities!=',num2str(max(abs(imag(nUpCal))))])
    end
    nUpCal = real(nUpCal);
    nDownCal = diag(ExpMatNorDownDownCal);      %%%%%%TESTING aug2021
    if max(abs(imag(nDownCal))) > 10^-15        
        disp(['imaginary densities!=',num2str(max(abs(imag(nDownCal))))])
    end
    nDownCal = real(nDownCal);    
% % % %     nUpCal = (abs(eVector(1:nBands,:)).^2)*fermi;%??????????????????????????
% % % %     nDownCal = (abs(eVector((3*nBands + 1):4*nBands,:)).^2)*(1-fermi);%???????????????????
    %anomalousncalc
    % problem in matlab R2018b and later
    if nargout_tmp > 6 
        szbdg=size(BdGMatrix);
    TotKE = 0.5*sum(diag((eVector)'*([BdGMatrix(1:szbdg(1)/2,1:szbdg(1)/2),zeros(szbdg(1)/2);...
        zeros(szbdg(1)/2),-BdGMatrix(szbdg(1)/2+1:szbdg(1),szbdg(1)/2+1:szbdg(1))])*(eVector)).*fermi);%.*tanh(En/kT));
    end
    % Mainak
%     uCal = (abs(eVector(1:nBands,floor(nBands/2)).^2));%*ones(size(fermi));
%     vCal = (abs(eVector((nBands + 1):end,floor(nBands/2))).^2);%*ones(size(1-fermi));
    % Mainak
