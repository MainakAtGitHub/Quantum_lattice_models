function KXX_summed=supcl_supfl(inputfile)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
read_input_file=inputfile;
read_input;
read_input_file


[filepath,name,ext]=fileparts(BdGfileName);
[~,nopath_casestring,caseext]=fileparts(casestring);
[~,nopath_Vimp,extVimp]=fileparts(Vimp);
LDOSfileName0 = [nopath_casestring,caseext,'_Vimp_',num2str(nopath_Vimp),num2str(extVimp),'_N_', num2str(N)];
if isempty(filepath)
    filepath='.';
end
dirstring=[filepath,'/','data_',LDOSfileName0];

kx_ind=[0:(M-1);ones(1,M)*M];
kx_ind=kx_ind./repmat(gcd(kx_ind(1,:),kx_ind(2,:)),2,1);
ky_ind=kx_ind;

load(TB_file);

if ~exist('pen_dep_dir','var')
    pen_dep_dir = 0;
end
if pen_dep_dir ~= 0
    [pen_dep_gamma,~] = make_pen_dep_gamma_supcl(N, TBparameters, latticeVector, pen_dep_dir);
    [pen_dep_gamma_dia,superLatticeVectors] = make_pen_dep_gamma_dia_supcl(N, TBparameters, latticeVector, pen_dep_dir);
end

kx = (2*pi/M)*(0:(M - 1));
ky = kx;

nSuperCells = size(superLatticeVectors,1);

KXX_summed=[0,0,0];
for index=1:M^2
    iKy= mod(index-1,M)+1;
    iKx= ceil(index/M);
    if iKy==1
        disp(['reading k-point',num2str(iKx),' ',num2str(iKy)])
    end;
    
    disp([iKx iKy]);
    k = [kx(iKx) ky(iKy)];
    kSpace_pen_dep_gamma = 0;
%     kSpace_pen_dep_gamma_c=0;
    kSpace_pen_dep_gamma_dia = 0;
%     kSpace_pen_dep_gamma_dia_c=0;
    
    for iUnitCell = 1:nSuperCells
        iLatticeVector = superLatticeVectors(iUnitCell,:);
        kSpace_pen_dep_gamma = kSpace_pen_dep_gamma + pen_dep_gamma(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
%         kSpace_pen_dep_gamma_c = kSpace_pen_dep_gamma_c + pen_dep_gamma(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
        kSpace_pen_dep_gamma_dia = kSpace_pen_dep_gamma_dia + pen_dep_gamma_dia(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
%         kSpace_pen_dep_gamma_dia_c = kSpace_pen_dep_gamma_dia_c + pen_dep_gamma_dia(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
    end

    ekukvk_fileGF=[dirstring,'/','kx_',num2str(kx_ind(1,iKx)),'_',num2str(kx_ind(2,iKx)),'ky_',num2str(ky_ind(1,iKy)),'_',num2str(ky_ind(2,iKy)),'_g_GF.mat'];
    if calcGreens
        load(ekukvk_fileGF,'Ek_vector','eigVectorK','eigValueK');
        
        En = eigValueK;
        eVector = eigVectorK;
        fermi = 1./(1 + exp(En/kT));
        
        DenMatrChk=(repmat(En,1,length(En)) - repmat(transpose(En),length(En),1));
        DenMatrSml=abs(DenMatrChk)<10^(-10);
        FracMatr= (repmat(fermi,1,length(fermi)) - repmat(transpose(fermi),length(fermi),1))./(repmat(En,1,length(En)) - repmat(transpose(En),length(En),1));
        RpMatr =  repmat(En,1,length(En));
        aux_mat = -1/kT*exp(RpMatr(DenMatrSml)/kT)./(1+exp(RpMatr(DenMatrSml)/kT)).^2;
        DenMatrBig = RpMatr(DenMatrSml)/kT > 50;
        aux_mat(DenMatrBig) = 0;
        FracMatr(DenMatrSml) = aux_mat;        
        
        KXX=sum(sum((abs(eVector'*[kSpace_pen_dep_gamma, zeros(size(kSpace_pen_dep_gamma)); zeros(size(kSpace_pen_dep_gamma)), kSpace_pen_dep_gamma]*eVector)).^2.*FracMatr))/(length(En)/2/nOrbitals);
        KXX_dia = sum(diag(eVector'*[kSpace_pen_dep_gamma_dia,zeros(size(kSpace_pen_dep_gamma_dia));zeros(size(kSpace_pen_dep_gamma_dia)),-kSpace_pen_dep_gamma_dia]*eVector).*fermi)/(length(En)/2/nOrbitals);
        
        KXX_tot = KXX_dia + KXX;
        KXX_summed=KXX_summed + [KXX, KXX_dia, KXX_tot];        
    end
end
KXX_summed=KXX_summed/M^2;
save([dirstring,'/',name,'_pen_dep_supcl_','M_',num2str(M),'_pen_dep_dir_',num2str(pen_dep_dir),'.txt'],'KXX_summed','-ascii');


