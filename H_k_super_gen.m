% H10,H20,.. = supercell_hopping;
% H=sum H10,.....fourier transformed;
% %     H = H +Himp;
% H = H+Hub_hund copied part;
% H = H+ Soc part;
% Pf([H(k), sum of 9 fourier transformed supercell_delta; (sum of 9 fourier transformed supercell_delta)', -transpose(H(-k))\tau_x])
%     ]
% H=H0+get_Himp;
% H = H+ hubbard_hund copied stff;
% H = H+ Soc terms;
% do H(k);

function pf_or_eig_val_list = H_k_super_gen(inputfile,klist, calc_pf,klist_no_points)

if ~isnumeric(calc_pf)
    calc_pf=str2num(calc_pf);
end
if nargin > 3
   if isstr(klist_no_points)
       klist_no_points=str2num(klist_no_points);
   end
end
if nargin < 4
    klist_no_points=100;
end
var_pos_neg = num2str(klist_no_points);
if ~isnumeric(klist)
    var_pos_neg=[var_pos_neg,klist];
    switch klist
        case 'y' 
            klist=[0*(0:1/klist_no_points:1);pi*(0:1/klist_no_points:1)]';
        case 'minus_y'
            klist=-[0*(0:1/klist_no_points:1);pi*(0:1/klist_no_points:1)]';
        case 'xy'
            klist=[pi*(0:1/klist_no_points:1);pi*(0:1/klist_no_points:1)]';
        otherwise
            disp('Not implemented');
    end
end

if nargin < 3
    calc_pf=false;
end
if calc_pf == 1
    klist=[0,0;...
              0,pi];
end

if calc_pf == 2
    klist=[0,0;...
              pi,pi];
end

     read_input_file=inputfile;
     read_input;
     read_input_file
load(TB_file);
mu=0;
load(BdGfileName);
maxHop = max(max(abs(latticeVector)));
nBands = N(1)*N(2)*nOrbitals;
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) = ...
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) - mu*eye(nOrbitals);

[H0_sup_hop, superLatticeVectors] = supercell_hoppings(N, TBparameters, latticeVector);
[deltaSuper,superDeltaVectors] = supercell_delta(nOrbitals, delta, maxHop,spin_and_nambu,N);
nSuperCells = size(superLatticeVectors,1);
if ~exist('klist','var')
    kx = 0;%(2*pi/M)*(0:(M - 1));
    ky = 0;%pi;
    
    klist = [kx,ky];%[kx(iKx) ky(iKy)];
end

if calc_pf
     pf_or_eig_val_list = zeros(length(klist(:,1)),1);
else
     pf_or_eig_val_list = zeros(N(1)*N(2)*nOrbitals*(1+spin_and_nambu)*2,length(klist(:,1)));
end
    
for itr_k =1:length(klist(:,1))
    k = klist(itr_k,:);
    kSpaceHopping = 0;
    kSpaceHoppingc=0;
    kSpaceGap = 0;
    for iUnitCell = 1:nSuperCells
        iLatticeVector = superLatticeVectors(iUnitCell,:);
        kSpaceHopping = kSpaceHopping + H0_sup_hop(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
        
        kSpaceHoppingc = kSpaceHoppingc + H0_sup_hop(:,:,iUnitCell)*exp(-1i*(iLatticeVector*k'));
        
        kSpaceGap = kSpaceGap + deltaSuper(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
    end
% H0_summed_sup_hop_k =  

HImpurity=get_Himp(Vimp,N,nOrbitals,sublattice,[],BdGfileName);

KESuper = kSpaceHopping + HImpurity;
KESuperdown = kSpaceHopping + HImpurity;
KESuperc = kSpaceHoppingc + HImpurity;
KESupercdown = kSpaceHoppingc + HImpurity;

KESuper=KESuper-field(end)*eye(nBands);
KESuperc=KESuperc-field(end)*eye(nBands);

KESuperdown=KESuperdown+field(end)*eye(nBands);
KESupercdown=KESupercdown+field(end)*eye(nBands);
% Hubbard-Hund and SOC below

KESuper = KESuper + U*diag(nDown);
KESuperc = KESuperc + U*diag(nDown);

for it1 = 1:length(nDown)
    IntUp1(it1) = U_pr*(sum(nDown([nOrbitals*(ceil(it1/nOrbitals)-1)+1:nOrbitals*ceil(it1/nOrbitals)]))-nDown(it1)) + (U_pr-J)*(sum(nUp([nOrbitals*(ceil(it1/nOrbitals)-1)+1:nOrbitals*ceil(it1/nOrbitals)]))-nUp(it1));
end
KESuper = KESuper + diag(IntUp1);
KESuperc = KESuperc + diag(IntUp1);

KESuperdown = KESuperdown + U*diag(nUp);
KESupercdown = KESupercdown + U*diag(nUp);
for it2 = 1:length(nUp)
    IntDown1(it2) = U_pr*(sum(nUp([nOrbitals*(ceil(it2/nOrbitals)-1)+1:nOrbitals*ceil(it2/nOrbitals)]))-nUp(it2)) + (U_pr-J)*(sum(nDown([nOrbitals*(ceil(it2/nOrbitals)-1)+1:nOrbitals*ceil(it2/nOrbitals)]))-nDown(it2));
end
KESuperdown = KESuperdown + diag(IntDown1);
KESupercdown = KESupercdown + diag(IntDown1);

if or(spinfullnormal,spin_and_nambu)
    H_off_up = -U*diag(nAnoUpDown);  %%23Mar2021, nAnoUpDown is expectation of cDownDaggercUp....**Note the flip of Up and Down, similar flip in nAnoUpDown
    H_off_down = -U*diag(nAnoDownUp); %%%??????????????COULD BE SETA AS conj(H_off_up)?
    if exist('field','var')
        if numel(field)>1
            H_off_up=H_off_up-field(1)*eye(nBands);
            H_off_down=H_off_down-field(1)*eye(nBands);
            
            H_off_up=H_off_up-(-1i)*field(2)*eye(nBands);
            H_off_down=H_off_down-(1i)*field(2)*eye(nBands);
        end
    end
    
    if int_soc     % SOC corresponding to the bilayer case in spinfullnormal........SOC corresponding to spin_and_nambu is encoded in variable "spin_and_nambu_SOC"
        H_off_up=H_off_up+H_soc1+(-1i)*H_soc2;
        H_off_down=H_off_down+H_soc1+(1i)*H_soc2;
    end
    if spin_and_nambu_SOC %% xz, yz, xy, x2-y2,z2 basis (%%%%%probably calculation of HSOC_UU etc can be done outside loop and then added here inside the loop, comment in May2021)
        lmbda_SOC=spin_and_nambu_SOC;%0.005;%0.01;%0.001;%0.02;
        %                         Lx5 = [[0 0 1i 0 0];[0 0 0 -1i -sqrt(3)*1i];[-1i 0 0 0 0];[0 1i 0 0 0];[0 sqrt(3)*1i 0 0 0]];
        %                         Ly5 = [[0 0 0 1i -sqrt(3)*1i];[0 0 1i 0 0];[0 -1i 0 0 0];[-1i 0 0 0 0];[sqrt(3)*1i 0 0 0 0]];
        %                         Lz5= [[0 -1i 0 0 0];[1i 0 0 0 0];[0 0 0 2i 0];[0 0 -2i 0 0];[0 0 0 0 0]];
        if ~exist('Orb_seq','var')
            % z^2 (1), xz (2), yz (3), x^2-y^2 (4), xy (5)
            % Orb_seq=[2,3,1,5,4];
            %Orb_seq=[2,3,4,5,1]; % correct seq for Ikeda for
            %our .csv file
            %Orb_seq=[3,2,1,5,4];
            %Orb_seq=[3,2,4,5,1];
            disp('Define orbital sequence for SO coupling')
        end
        Lx5 = [[0             0  sqrt(3)*1i      0    0];...
            [0             0           0      0   1i];...
            [-sqrt(3)*1i   0           0    -1i    0];...
            [0             0          1i      0    0];...
            [0           -1i           0      0    0]];
        
        Ly5 = [[0    sqrt(3)*1i           0      0    0];...
            [-sqrt(3)*1i   0           0     1i    0];...
            [0             0           0      0   1i];...
            [0           -1i           0      0    0];...
            [0             0         -1i      0    0]];
        
        Lz5 = [[0             0           0      0    0];...
            [0             0         -1i      0    0];...
            [0            1i           0      0    0];...
            [0             0           0      0  -2i];...
            [0             0           0     2i    0]];
        
        Lx5=Lx5(:,Orb_seq); Lx5= Lx5(Orb_seq,:);
        Ly5=Ly5(:,Orb_seq); Ly5= Ly5(Orb_seq,:);
        Lz5=Lz5(:,Orb_seq); Lz5= Lz5(Orb_seq,:);
        
        HSOC_UU = Lz5(1:nOrbitals,1:nOrbitals)/2;
        HSOC_DD = -Lz5(1:nOrbitals,1:nOrbitals)/2;
        HSOC_UD = (Lx5(1:nOrbitals,1:nOrbitals)+1i*Ly5(1:nOrbitals,1:nOrbitals))/2;
        HSOC_DU = (Lx5(1:nOrbitals,1:nOrbitals)-1i*Ly5(1:nOrbitals,1:nOrbitals))/2;
        
        H_off_up=H_off_up+lmbda_SOC*kron(eye(N(1)*N(2)),HSOC_UD);%repmat(HSOC_UD,N^2);
        H_off_down=H_off_down+lmbda_SOC*kron(eye(N(1)*N(2)),HSOC_DU);%repmat(HSOC_DU,N^2);
        KESuper = KESuper + lmbda_SOC*kron(eye(N(1)*N(2)),HSOC_UU);%repmat(HSOC_UU,N^2);
        KESuperc = KESuperc + lmbda_SOC*kron(eye(N(1)*N(2)),HSOC_UU);%repmat(HSOC_UU,N^2);
        KESuperdown = KESuperdown + lmbda_SOC*kron(eye(N(1)*N(2)),HSOC_DD);%repmat(HSOC_DD,N^2);
        KESupercdown = KESupercdown + lmbda_SOC*kron(eye(N(1)*N(2)),HSOC_DD);%repmat(HSOC_DD,N^2);
    end
    E_Hund_remaining=0; % initializing E_Hund_remaining
    for lat_pt_num = 1:nBands/nOrbitals % form the Hund interaction using the ExpMatNor matrix to be taken as an output (just like nUp, nDown etc.) from BdG_step2 fn
        if ~exist('J_pr','var') %temporary to avoid error, can be removed later
            J_pr=0;
        end
        
        for orb_num1 = 1:nOrbitals
            
            dmi=(lat_pt_num-1)*nOrbitals+orb_num1; %calculating matrix indices
            
            temp_diag=diag(ExpMatNorDownUp);
            H_Hund_UpDown(dmi,dmi) = ...; % forming upper diagonal part of the matrix H_Hund_UpDown (only J because U has been implemented in earlier step)
                (-J)*(sum(temp_diag((lat_pt_num-1)*nOrbitals+1:(lat_pt_num-1)*nOrbitals+nOrbitals))-temp_diag(dmi));
            
            E_Hund_remaining=E_Hund_remaining + H_Hund_UpDown(dmi,dmi)*ExpMatNorUpDown(dmi,dmi); %%%%%% FIRST OCCURENCE OF E_Hund_remaining, SO NO "E_Hund_remaining=E_Hund_remaining+"
            
            temp_diag=diag(ExpMatNorUpDown);
            H_Hund_DownUp(dmi,dmi) = ...; % forming upper diagonal part of the matrix H_Hund_DownUp (only J because U has been implemented in earlier step)
                (-J)*(sum(temp_diag((lat_pt_num-1)*nOrbitals+1:(lat_pt_num-1)*nOrbitals+nOrbitals))-temp_diag(dmi));
            
            E_Hund_remaining=E_Hund_remaining + H_Hund_DownUp(dmi,dmi)*ExpMatNorDownUp(dmi,dmi);
            
            if orb_num1 < nOrbitals
                for orb_num2 = orb_num1+1:nOrbitals
                    
                    dmj=(lat_pt_num-1)*nOrbitals+orb_num2; %calculating matrix indices
                    
                    H_Hund_UpUp(dmi,dmj) = ... % forming upper triangular part of the matrix H_Hund_UpUp
                        J_pr*ExpMatNorDownDown(dmi,dmj)...
                        -(U_pr-J)*ExpMatNorUpUp(dmj,dmi)...
                        +J*ExpMatNorDownDown(dmj,dmi);
                    
                    E_Hund_remaining=E_Hund_remaining + H_Hund_UpUp(dmi,dmj)*ExpMatNorUpUp(dmi,dmj);
                    
                    H_Hund_UpUp(dmj,dmi) = ... % forming lower triangular part of the matrix H_Hund_UpUp
                        J_pr*ExpMatNorDownDown(dmj,dmi)...
                        -(U_pr-J)*ExpMatNorUpUp(dmi,dmj)...
                        +J*ExpMatNorDownDown(dmi,dmj);
                    
                    E_Hund_remaining=E_Hund_remaining + H_Hund_UpUp(dmj,dmi)*ExpMatNorUpUp(dmj,dmi);
                    
                    H_Hund_DownDown(dmi,dmj) = ... % forming upper triangular part of the matrix H_Hund_DownDown
                        J_pr*ExpMatNorUpUp(dmi,dmj)...
                        -(U_pr-J)*ExpMatNorDownDown(dmj,dmi)...
                        +J*ExpMatNorUpUp(dmj,dmi);
                    
                    E_Hund_remaining=E_Hund_remaining + H_Hund_DownDown(dmi,dmj)*ExpMatNorDownDown(dmi,dmj);
                    
                    H_Hund_DownDown(dmj,dmi) = ... % forming lower triangular part of the matrix H_Hund_DownDown
                        J_pr*ExpMatNorUpUp(dmj,dmi)...
                        -(U_pr-J)*ExpMatNorDownDown(dmi,dmj)...
                        +J*ExpMatNorUpUp(dmi,dmj);
                    
                    E_Hund_remaining=E_Hund_remaining + H_Hund_DownDown(dmj,dmi)*ExpMatNorDownDown(dmj,dmi);
                    
                    H_Hund_UpDown(dmi,dmj) = ... % forming upper triangular part of the matrix H_Hund_UpDown
                        (-U_pr)*ExpMatNorDownUp(dmj,dmi)...
                        -J_pr*ExpMatNorDownUp(dmi,dmj);
                    
                    E_Hund_remaining=E_Hund_remaining + H_Hund_UpDown(dmi,dmj)*ExpMatNorUpDown(dmi,dmj);
                    
                    H_Hund_UpDown(dmj,dmi) = ... % forming lower triangular part of the matrix H_Hund_UpDown
                        (-U_pr)*ExpMatNorDownUp(dmi,dmj)...
                        -J_pr*ExpMatNorDownUp(dmj,dmi);
                    
                    E_Hund_remaining=E_Hund_remaining + H_Hund_UpDown(dmj,dmi)*ExpMatNorUpDown(dmj,dmi);
                    
                    H_Hund_DownUp(dmi,dmj) = ... % forming upper triangular part of the matrix H_Hund_DownUp
                        (-U_pr)*ExpMatNorUpDown(dmj,dmi)...
                        -J_pr*ExpMatNorUpDown(dmi,dmj);
                    
                    E_Hund_remaining=E_Hund_remaining + H_Hund_DownUp(dmi,dmj)*ExpMatNorDownUp(dmi,dmj);
                    
                    H_Hund_DownUp(dmj,dmi) = ... % forming lower triangular part of the matrix H_Hund_DownUp
                        (-U_pr)*ExpMatNorUpDown(dmi,dmj)...
                        -J_pr*ExpMatNorUpDown(dmj,dmi);
                    
                    E_Hund_remaining=E_Hund_remaining + H_Hund_DownUp(dmj,dmi)*ExpMatNorDownUp(dmj,dmi);
                end
            end
        end
    end
    
    KESuper = KESuper + H_Hund_UpUp;
    KESuperc = KESuperc + H_Hund_UpUp;
    KESuperdown = KESuperdown + H_Hund_DownDown;
    KESupercdown = KESupercdown + H_Hund_DownDown;
    H_off_up = H_off_up + H_Hund_UpDown;
    H_off_down = H_off_down + H_Hund_DownUp;
    
end

kSpaceHamiltonian = [[KESuper,H_off_up;H_off_down,KESuperdown], -kSpaceGap; -kSpaceGap', -conj([KESuperc,H_off_up;H_off_down,KESupercdown])];
testvr=max(max(abs(kSpaceHamiltonian-kSpaceHamiltonian')))
kSpaceHamiltonian=(kSpaceHamiltonian+kSpaceHamiltonian')/2;



if calc_pf
%     kSpaceHamiltonian=kSpaceHamiltonian./(abs(prod(diag(evls)))).^(1/numel(diag(evls)));
    H_Pf_k = kSpaceHamiltonian*[zeros(size(kSpaceGap)), eye(size(kSpaceGap)); eye(size(kSpaceGap)),zeros(size(kSpaceGap))];
    pf_or_eig_val_list(itr_k)=pfaffian_LTL_AK_MP(H_Pf_k);
else
    [evcs,evls]= eig(kSpaceHamiltonian);
    pf_or_eig_val_list(:,itr_k)=(diag(evls))';
end
disp(['itr_k=',num2str(itr_k)]);
% top_inv=pfaffian_LTL(H_Pf_k);
end
[pthstr,nm,xt] = fileparts(inputfile); 

if ~isempty(pthstr)
    pthstr=[pthstr,'/'];
end
if calc_pf
    save([pthstr,'pfaffian_sign_at_0_pi.mat'],'pf_or_eig_val_list');
    save([pthstr,'pfaffian_sign_at_0_pi.txt'],'pf_or_eig_val_list','-ascii');
else
    save([pthstr,'sup_cell_k_spectr_',var_pos_neg,'.mat'],'klist','pf_or_eig_val_list');
%   save([pthstr,'/pfaffian_sign_at_0_pi.txt'],'klist','pf_or_eig_val_list','-ascii');
end
% if calc_pf
%     topo_inv = prod(pf_or_eig_val_list);
% else
%     topo_inv;
% end


