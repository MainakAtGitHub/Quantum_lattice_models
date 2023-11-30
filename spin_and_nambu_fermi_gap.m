function [ output_args ] = spin_and_nambu_fermi_gap(inputfile,k_points,onsite_SOC,switch_on_nm_dg_br,do_rot_d_vec)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if nargin<5
    do_rot_d_vec=false;
end

if nargin <4
    switch_on_nm_dg_br=0; % switch on numerical degeneracy breaker
end
if nargin < 3
    onsite_SOC=0;  % onsite_SOC used as boolean
end

if ~exist('spin_and_nambu','var')
    spin_and_nambu=0;
end
read_input_file=inputfile;
read_input;
read_input_file

load(TB_file,'-mat');
nOrbitals = size(TBparameters,1);
load(Gamma_file,'-mat');
load(BdGfileName,'-mat');
% N = sqrt(size(delta,1)/(nOrbitals*(spin_and_nambu+1)));
if numel(N)==1
    N=[N,N];
end
nUnitCellsDelta = size(latticeVectorsSC,1);

deltaCenter = zeros(nOrbitals*(spin_and_nambu+1), nOrbitals*(spin_and_nambu+1), nUnitCellsDelta);
jCell = [ceil(N(1)/2) ceil(N(2)/2)];
for i = 1:nUnitCellsDelta
    iCell = jCell + latticeVectorsSC(i,:);
    [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, iCell, jCell);
    if ~(spin_and_nambu)
        deltaCenter(:,:,i) = delta(iRange, jRange);
    else
        dlta1= delta(1:size(delta,1)/2,1:size(delta,1)/2);
        dlta2= delta(1:size(delta,1)/2,size(delta,1)/2+1:size(delta,1));
        dlta3= delta(size(delta,1)/2+1:size(delta,1),1:size(delta,1)/2);
        dlta4= delta(size(delta,1)/2+1:size(delta,1),size(delta,1)/2+1:size(delta,1));
        deltaCenter(:,:,i) = [dlta1(iRange, jRange),dlta2(iRange, jRange);...
            dlta3(iRange, jRange),dlta4(iRange, jRange)];
    end
end
delta = deltaCenter;



% set up a mesh in the BZ
k_range=[-pi:pi/k_points:pi];
[kx,ky]=meshgrid(k_range);
% put the k-points into the correct vectorlike format
klist=[kx(:),ky(:)];

% calculate the eigenvalues in the normal state for all k points
% if ~onsite_SOC
% [kSpaceEigenValuesNormal,kSpaceEigenVectorsNormal]=bands_precalc_klist2(inputfile,klist);
% end
% if onsite_SOC
    [kSpaceEigenValuesNormal,kSpaceEigenVectorsNormal]=bands_precalc_klist2(inputfile,klist,onsite_SOC); %%###Don't consider SOC for extracting Fermi surface? %yes, do consider soc
% end



whole_fermi_contour=[];
band_no = [];
no_of_bands=nOrbitals*(onsite_SOC+1);

% colorlist = [1,0,0; 0.6,0,0; 0,1,0; 0,0.6,0; 0,0,1; 0,0, 0.6; 1,0,1; 0.6,0,0.6; 0,0.6,0.6; 0,0.6,0.6]; %r g b magenta cyan and lighter shades
colorlist = [1,0,0; 0,1,0; 0,0,1; 1,0,1; 0,1,1]; % rgbmc
% colorlist = [colorlist(1:2:end,:);colorlist(2:2:end,:)]; %% Just shuffling the colorlist to put the spin ups color together and the spin downs too for the same colors
% colorlist = [1,0,0; 0,1,0; 0,0,1; 1,0,1; 0,1,1];
% colorlist = [] % colors taken from "http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html#2"

% colorlist = ['r','g','b','c','m'];
weights=[];

for j=1:(1+onsite_SOC):no_of_bands
    eigv_reshaped = reshape(kSpaceEigenValuesNormal(:,j),length(k_range),length(k_range));
    C=contourc(k_range,k_range,eigv_reshaped,[-0,-0]);
    contour_coordinates{j}=[];
    level_index=1;
    while level_index < size(C,2)
        level_point_no=C(2,level_index);
        contour_coordinates{j}=[contour_coordinates{j},[C(1,level_index+1:level_index+level_point_no);C(2,level_index+1:level_index+level_point_no)]];
        level_index=level_index+level_point_no+1;
    end
    if ~isempty(contour_coordinates{j})
        whole_fermi_contour=[whole_fermi_contour,contour_coordinates{j}];
        band_no=[band_no,j*ones(1,size(contour_coordinates{j},2))];
    end
end



[kSpaceEigenValuesNormal_fermi,kSpaceEigenVectorsNormal_fermi]=bands_precalc_klist2(inputfile,whole_fermi_contour',onsite_SOC,switch_on_nm_dg_br);
[m_kSpaceEigenValuesNormal_fermi,m_kSpaceEigenVectorsNormal_fermi]=bands_precalc_klist2(inputfile,-whole_fermi_contour',onsite_SOC,switch_on_nm_dg_br);

if onsite_SOC
%#    kSpaceEigenValuesNormal_fermi = [kSpaceEigenValuesNormal_fermi(:,1:2:end),kSpaceEigenValuesNormal_fermi(:,2:2:end)]; %%%%%% putting up spin eigvals and
%#    kSpaceEigenVectorsNormal_fermi = kSpaceEigenVectorsNormal_fermi(:,:,[1:2:no_of_bands,2:2:no_of_bands]);
%     kSpaceEigenVectorsNormal_fermi(:,:,1:end/2) = kSpaceEigenVectorsNormal_fermi_dummy(:,:,1:2:end);  %%%%% eigvecs together and down too
%     kSpaceEigenVectorsNormal_fermi(:,:,end/2+1:end) = kSpaceEigenVectorsNormal_fermi_dummy(:,:,2:2:end);
    
%#    m_kSpaceEigenValuesNormal_fermi = [m_kSpaceEigenValuesNormal_fermi(:,1:2:end),m_kSpaceEigenValuesNormal_fermi(:,2:2:end)];   %%%%%% putting up spin eigvals and
%#    m_kSpaceEigenVectorsNormal_fermi = m_kSpaceEigenVectorsNormal_fermi(:,:,[1:2:no_of_bands,2:2:no_of_bands]);
%     m_kSpaceEigenVectorsNormal_fermi(:,:,1:end/2) = m_kSpaceEigenVectorsNormal_fermi_dummy(:,:,1:2:end); %%%%% eigvecs together and down too
%     m_kSpaceEigenVectorsNormal_fermi(:,:,end/2+1:end) = m_kSpaceEigenVectorsNormal_fermi_dummy(:,:,2:2:end);
    
    clear kSpaceEigenVectorsNormal_fermi_dummy m_kSpaceEigenVectorsNormal_fermi_dummy;
    
    colorlist = [colorlist;colorlist]; 
    
%     for i = 1:2:no_of_bands % shuffling the band numbers too for putting spin up bands together and spin down bands together
%         new_band_no(find(band_no==i))=ceil(i/2);
%     end
    
%     for i = 2:2:no_of_bands % shuffling the band numbers too for putting spin up bands together and spin down bands together
%         new_band_no(find(band_no==i))=(no_of_bands+i)/2;
%     end
%     band_no = new_band_no; clear new_band_no;
end

kSpaceGap_memory=zeros(size(whole_fermi_contour',1),2*nOrbitals,2*nOrbitals);
for iK = 1:size(whole_fermi_contour',1)
    
    skz_kSpaceEigenVectorsNormal_fermi = squeeze(kSpaceEigenVectorsNormal_fermi(iK,:,:));
    skz_m_kSpaceEigenVectorsNormal_fermi = squeeze(m_kSpaceEigenVectorsNormal_fermi(iK,:,:));
    
    skz_m_kSpaceEigenVectorsNormal_fermi=conj(skz_m_kSpaceEigenVectorsNormal_fermi);
    
    vectr=abs(skz_kSpaceEigenVectorsNormal_fermi(band_no(iK),:));
    
    if onsite_SOC
    halved_squared_vectr=(vectr(1:no_of_bands/2)).^2+(vectr(no_of_bands/2+1:end)).^2;
    else
        halved_squared_vectr=vectr.^2;
    end
    
    orbtl_no = max(find(halved_squared_vectr==max(halved_squared_vectr)));
    weights = [weights;colorlist(ceil(orbtl_no/(1)),:)]; %+onsite_SOC
    
    kSpaceGap = zeros((spin_and_nambu+1)*nOrbitals,(spin_and_nambu+1)*nOrbitals); % (spin_and_nambu+1) for distinguishing 2x2 and 4x4 spin_and_nambu case
    
            for iUnitCellDelta = 1:nUnitCellsDelta
                iLatticeVectorDelta = latticeVectorsSC(iUnitCellDelta,:);
                kSpaceGap = kSpaceGap + delta(:,:,iUnitCellDelta)*exp(1i*([iLatticeVectorDelta*whole_fermi_contour(:,iK)]'));
            end
            
    if ~onsite_SOC 
    band_gap=[skz_kSpaceEigenVectorsNormal_fermi', zeros(size(skz_kSpaceEigenVectorsNormal_fermi));...
        zeros(size(skz_kSpaceEigenVectorsNormal_fermi)), skz_kSpaceEigenVectorsNormal_fermi']'...
        *kSpaceGap*...
        conj([skz_m_kSpaceEigenVectorsNormal_fermi',zeros(size(skz_m_kSpaceEigenVectorsNormal_fermi));...
        zeros(size(skz_m_kSpaceEigenVectorsNormal_fermi)),skz_m_kSpaceEigenVectorsNormal_fermi']);
    end 
    if onsite_SOC 
       skz_kSpaceEigenVectorsNormal_fermi=skz_kSpaceEigenVectorsNormal_fermi([1:2:2*nOrbitals,2:2:2*nOrbitals],:); 
       skz_m_kSpaceEigenVectorsNormal_fermi=skz_m_kSpaceEigenVectorsNormal_fermi([1:2:2*nOrbitals,2:2:2*nOrbitals],:);
       band_gap=skz_kSpaceEigenVectorsNormal_fermi*kSpaceGap*skz_m_kSpaceEigenVectorsNormal_fermi';  
    end
    kSpaceGap_memory(iK,:,:)=band_gap;
%     proj_band_gap = band_gap([ceil(band_no(iK)/(1)),ceil(band_no(iK)/(1))+no_of_bands/(1)],[ceil(band_no(iK)/(1)),ceil(band_no(iK)/(1))+no_of_bands/(1)]);
%       if mod(orbtl_no,5) ~= 0
%           orbtl_no_mod_5 = mod(orbtl_no,5);
%       else 
%           orbtl_no_mod_5=5;
%       end
%       orbtl_no=orbtl_no_mod_5;
      
%       proj_band_gap = band_gap([ceil(band_no(iK)/(1+onsite_SOC)),ceil(band_no(iK)/(1+onsite_SOC))+no_of_bands/(1+onsite_SOC)],...
%           [ceil(band_no(iK)/(1+onsite_SOC)),ceil(band_no(iK)/(1+onsite_SOC))+no_of_bands/(1+onsite_SOC)]);
      if onsite_SOC
      band_no_iK=(band_no(iK)+1)/2;    
      proj_band_gap = band_gap([band_no_iK,band_no_iK+nOrbitals],[band_no_iK,band_no_iK+nOrbitals]);
      else
          proj_band_gap = band_gap([band_no(iK),band_no(iK)+no_of_bands],[band_no(iK),band_no(iK)+no_of_bands]);
      end
    
    d0 = (proj_band_gap(1,2) - proj_band_gap(2,1))/2;
    d3 = (proj_band_gap(1,2) + proj_band_gap(2,1))/2;
    d1 = (proj_band_gap(2,2) - proj_band_gap(1,1))/2;
    d2 = (proj_band_gap(2,2) + proj_band_gap(1,1))/2/1i;
    
    %%%%%%%%% transformation to d vectors;
      d_vecs(iK,:)=[d0,d1,d2,d3];
%     d_vecs(iK,:)=[d0,d1,d2,d3];
end

        load(BdGfileName);
        [filepath,name,ext]=fileparts(BdGfileName);
       

figure;
% stem3(whole_fermi_contour(1,:),whole_fermi_contour(2,:),d_vecs(:,1),'MarkerFaceColor','r');
scatter(whole_fermi_contour(1,:)/pi,whole_fermi_contour(2,:)/pi,[],real(d_vecs(:,1)),'.');
blue_red_map(gcf);
axis square;
ca=caxis;
mx_ca=max(abs(ca));
caxis([-mx_ca,mx_ca]);
colorbar;
xlabel('k_x/pi');ylabel('k_y/pi'),zlabel('singlet band gap');
title(['singlet bandgap (real)']);
testgca=gca;
testgca.FontSize=15;
box on;
axis square;
xlim([-1,1]);
ylim([-1,1]);
 print_pdf([filepath,filesep,'singlet bandgap (real)']);
 saveas(gcf,[filepath,filesep,'singlet bandgap (real)']);

figure;
% stem3(whole_fermi_contour(1,:),whole_fermi_contour(2,:),d_vecs(:,1),'MarkerFaceColor','r');
scatter(whole_fermi_contour(1,:)/pi,whole_fermi_contour(2,:)/pi,[],imag(d_vecs(:,1)),'.');
blue_red_map(gcf);
axis square;
ca=caxis;
mx_ca=max(abs(ca));
caxis([-mx_ca,mx_ca]);
colorbar;
xlabel('k_x/pi');ylabel('k_y/pi'),zlabel('singlet band gap');
title(['singlet bandgap (imag)']);
testgca=gca;
testgca.FontSize=15;
box on;
axis square;
xlim([-1,1]);
ylim([-1,1]);
%hold on;
 print_pdf([filepath,filesep,'singlet bandgap (imag)']);
  saveas(gcf,[filepath,filesep,'singlet bandgap (imag)']);

figure;
% stem3(whole_fermi_contour(1,:),whole_fermi_contour(2,:),d_vecs(:,1),'MarkerFaceColor','r');
mag_d_vec=sqrt(sum(d_vecs.*conj(d_vecs),2));
scatter(whole_fermi_contour(1,:)/pi,whole_fermi_contour(2,:)/pi,[],mag_d_vec,'.');
blue_red_map(gcf);
axis square;
ca=caxis;
mx_ca=max(abs(ca));
caxis([-mx_ca,mx_ca]);
colorbar;
xlabel('k_x/pi');ylabel('k_y/pi'),zlabel('absolute band gap');
title(['absolute bandgap']);
testgca=gca;
testgca.FontSize=15;
box on;
axis square;
xlim([-1,1]);
ylim([-1,1]);
 print_pdf([filepath,filesep,'absolute bandgap']);
 saveas(gcf,[filepath,filesep,'absolute bandgap']);

d11=(d_vecs(1,1));
d12=(d_vecs(1,2));
d13=(d_vecs(1,3));
d14=(d_vecs(1,4));
rot_d_vec=[d_vecs(:,1)/exp(1i*atan2(imag(d11),real(d11))), d_vecs(:,2)/exp(1i*atan2(imag(d12),real(d12))), d_vecs(:,3)/exp(1i*atan2(imag(d13),real(d13))) , d_vecs(:,4)/exp(1i*atan2(imag(d14),real(d14)))];
if do_rot_d_vec==true
    d_vecs=rot_d_vec;
end


figure;
x=whole_fermi_contour(1,:)/pi;
y=whole_fermi_contour(2,:)/pi;
z=0*whole_fermi_contour(1,:);
u=d_vecs(:,2);
v=d_vecs(:,3);
w=d_vecs(:,4);
q=quiver3(x',y',z',real(u),real(v),real(w));
hold on;
plot3(x',y',z','.');
        camlight head;
        lighting phong;
% c=q.Color;
% q.Color = band_no;
xlabel('k_x/pi');ylabel('k_y/pi'),zlabel('triplet band gap (real)');
% hold on;
% surf(kx,ky,zeros(size(kx)));
title({[' real u_{max}=',num2str(max(real(u))),' real v_{max}=',num2str(max(real(v))),' real w_{max}=',num2str(max(real(w)))]});
testgca=gca;
testgca.FontSize=15;
box on;
axis square;
xlim([-1,1]);
ylim([-1,1]);
view(2)
 print_pdf([filepath,filesep,'triplet band gap (real)']);
 saveas(gcf,[filepath,filesep,'triplet band gap (real)']);

figure;
% x=whole_fermi_contour(1,:);
% y=whole_fermi_contour(2,:);
% z=0*whole_fermi_contour(1,:);
% u=d_vecs(:,2);
% v=d_vecs(:,3);
% w=d_vecs(:,4);
q=quiver3(x',y',z',imag(u),imag(v),imag(w));
hold on;
plot3(x',y',z','.');
        camlight head;
        lighting phong;
% c=q.Color;
% q.Color = band_no;
xlabel('k_x/pi');ylabel('k_y/pi'),zlabel('triplet band gap (imag)');
% hold on;
% surf(kx,ky,zeros(size(kx)));
title({[' imag u_{max}=',num2str(max(imag(u))),' imag v_{max}=',num2str(max(imag(v))),' imag w_{max}=',num2str(max(imag(w)))]});
testgca=gca;
testgca.FontSize=15;
box on;
axis square;
xlim([-1,1]);
ylim([-1,1]);
view(2)
 print_pdf([filepath,filesep,'triplet band gap (imag)']);
 saveas(gcf,[filepath,filesep,'triplet band gap (imag)']);


figure;
scatter(whole_fermi_contour(1,:)/pi,whole_fermi_contour(2,:)/pi,[],weights,'.');
% blue_red_map(gcf);
axis square;
% ca=caxis;
% mx_ca=max(abs(ca));
% caxis([-mx_ca,mx_ca]);
% colorbar;
xlabel('k_x');ylabel('k_y'),zlabel('orbital weights');
title('orbital weights (colorbar)')
testgca=gca;
testgca.FontSize=15;
box on;
axis square;
xlim([-1,1]);
ylim([-1,1]);
 print_pdf([filepath,filesep,'orbital weights']);
 saveas(gcf,[filepath,filesep,'orbital weights']);
end

