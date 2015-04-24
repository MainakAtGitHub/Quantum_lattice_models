function p=process_Gamma(flnm,nOrbitals,sublattice)
% read Gamma from file
if nargin <2
    nOrbitals = 10;
end;
if nargin <3
    sublattice=1;
end;
nSites = 11;
list_plane_small_l=load('list_plane_small_l.csv');
if nargin < 1
    %flnm='ChiqData_LiFeAs_5_ARPES_v2_2Dpi.dat_rlist';
    %flnm = 'Gamma_fese_Toms_BS.dat_rlist';
    flnm = 'Gamma_Tom_U_0.95.dat_rlist';
end;
r=read_complex_matrix(flnm,nOrbitals^2,nOrbitals^2,2+1i);
r = real(r);

% extract Intra and Mixed pairing vertices (\Gamma_1111 & \Gamma_1221)
GammaMixed = zeros(nOrbitals^2,nSites^2);
count = 0;
for i = 1:nOrbitals
    for j = 1:nOrbitals
        count = count + 1;
        GammaMixed(count,:) = r((i-1)*nOrbitals^3 + (i-1)*nOrbitals^2 + (j-1)*nOrbitals + (j-1) + 1,:);
    end
end

latticeVectors = list_plane_small_l(:,1:2);
Gamma = reshape(GammaMixed, nOrbitals, nOrbitals, nSites^2);
orb2=nOrbitals/2;
if abs(sublattice)>0
    orb2=nOrbitals/2;
    GammaNewOrder = Gamma;
    GammaNewOrder(1:orb2,orb2+1:nOrbitals,:) = Gamma(orb2+1:nOrbitals,1:orb2,:);
    GammaNewOrder(orb2+1:nOrbitals,1:orb2,:) = Gamma(1:orb2,orb2+1:nOrbitals,:);
    Gamma = GammaNewOrder;
%elseif sublattice<0
%    orb2=nOrbitals/2;
else
    orb2=nOrbitals;
end;

% saving

latticeVectorsSC = latticeVectors;
Gamma = - Gamma; % change potential's sign convention
save([flnm,'.mat'], 'Gamma', 'latticeVectorsSC');

% % Shortening the range
cut = input('enter the Gamma cut   ');
latticeVectorsSCCut = zeros((2*cut + 1)^2,2);
GammaCut = zeros(nOrbitals, nOrbitals, (2*cut + 1)^2);
count = 0;
for i = -cut:cut
    for j = -cut:cut
        count = count + 1;
        latticeVectorsSCCut(count,:) = [i j];
        GammaCut(:,:,count) = Gamma(:,:,((latticeVectorsSC(:,1) == i) & (latticeVectorsSC(:,2) == j)));
        % special care for the edges, assume certain relation of Fe1 and
        % Fe2 in elementary cell.
        if abs(sublattice)>0
        if ((i==-cut) | (j==cut))
            GammaCut(1:nOrbitals/2,nOrbitals/2+1:nOrbitals,count)=zeros(nOrbitals/2);
        end;
        if ((i==cut) | (j==-cut))
            GammaCut(nOrbitals/2+1:nOrbitals,1:nOrbitals/2,count)=zeros(nOrbitals/2);
        end
        end;
    end
end
Gammatemp=Gamma;
latticeVectorsSCtemp=latticeVectorsSC;
Gamma=GammaCut;
latticeVectorsSC=latticeVectorsSCCut;
save([flnm,'cut',num2str(cut),'.mat'], 'Gamma', 'latticeVectorsSC');
Gamma=Gammatemp;
latticeVectorsSC=latticeVectorsSCtemp;

% plotting
%convert 2Fe to 1Fe
N = nSites;
if abs(sublattice) >0
centerCell = [ceil(N/2) ceil(N/2)];
latticeVectors1Fe = [];%zeros(2,2*N^2);
Gamma1Fe = [];%zeros(nOrbitals, nOrbitals, 2*N^2);
count = 0;
for i = -(ceil(N/2)-1):(ceil(N/2)-1)
    for j = -(ceil(N/2)-1):(ceil(N/2)-1)
        count = count + 1;
        Gamma2Fe = Gamma(:,:,count);
        switch sublattice
            case 1
                latticeVectors1Fe = [latticeVectors1Fe; i+j j-i; i+j j-i+1];     
            case -1
%                latticeVectors1Fe = [latticeVectors1Fe; i+j j-i+1; i+j j-i];
                        latticeVectors1Fe = [latticeVectors1Fe; i+j j-i+1; i+j j-i];
        end;  
        if sublattice >0
            Gamma1Fe = [Gamma1Fe Gamma2Fe(1:nOrbitals/2,:)];
        else
                      %  Gamma1Fe = [Gamma1Fe Gamma2Fe(1:nOrbitals/2,:)];
            Gamma1Fe = [Gamma1Fe Gamma2Fe(nOrbitals/2+1:nOrbitals,:)];
        end;

    end
end
Gamma1Fe =  reshape(Gamma1Fe,nOrbitals/2, nOrbitals/2, 2*N^2);
else
    Gamma1Fe=Gamma;
    latticeVectors1Fe=latticeVectorsSC;
end;
% remove on site potentials
 Gamma1Fe(:,:,(latticeVectors1Fe(:,1) == 0) & (latticeVectors1Fe(:,2) == 0))= 0;
Gamma2Plot = zeros(orb2*N, orb2*N);
for iOrbital = 1:orb2
    for jOrbital = 1:orb2
        GammaBlock = zeros(N,N);
        for i = 1:N
            for j = 1:N
                latticeVector = [j-ceil(N/2), ceil(N/2)-i];
                ind = find((latticeVectors1Fe(:,1) == latticeVector(1)) & (latticeVectors1Fe(:,2) == latticeVector(2)));
                GammaBlock(i,j) = Gamma1Fe(iOrbital, jOrbital, ind);
            end
        end
        Gamma2Plot(((iOrbital-1)*N + 1):iOrbital*N, ((jOrbital-1)*N + 1):jOrbital*N ) = GammaBlock;
    end
end
%numl = N;
tickx={'$d_{z^2}$','$d_{x^2-y^2}$','$d_{yz}$','$d_{xz}$','$d_{xy}$'}; % orbitals order
% global colorred
% colorred=false; % set colorscale
% figure1 = figure;
% global fsz
% fsz = 25; % font
% GammaRealMax = max(abs(real(Gamma2Plot(:)))); % for setting colormap
% scale = 's';
% if scale == 's'
%     image(sign(real(Gamma2Plot)).*sqrt(abs(real(Gamma2Plot)/GammaRealMax))*128+128);
% else
%     image(sign(real(Gamma2Plot)).*(abs(real(Gamma2Plot)/GammaRealMax))*128+128);
% end

r=realspaceplot(Gamma2Plot,N,tickx,flnm,'s');

% 
% scale='s';
% figure1=figure
% global colorred
% colorred=false; % set colorscale
% global fsz
% fsz = 25; % font
%  GammaRealMax = max(abs(real(Gamma2Plot(:)))); % for setting colormap
% % if scale == 's'
% %     image(sign(real(delta2Plot)).*sqrt(abs(real(delta2Plot)/deltaRealMax))*128+128);
% % else
% %     image(sign(real(delta2Plot)).*(abs(real(delta2Plot)/deltaRealMax))*128+128);
% % end
% % colorbar_rwb(figure1,deltaRealMax,ticks);
% % label_boxes(nOrbitals/2,numl,tickx);
% 
% % set the scale
% lm=log(GammaRealMax)/log(10);
% mtix=10^(ceil(lm));
% % do some refinement to avoid only single labels
% if (ceil(lm)-lm > 0.5)
%     tx=[-.5:.05:.5];
% elseif    (ceil(lm)-lm > 0.2)
%     tx=[-1:0.1:1];
% else
%     tx=[-1:0.2:1];
% end;
% ticks=mtix*tx; 
% labels = num2str(repmat(sign(ticks).*(abs(ticks)), 1, 1)', 2);
% if scale=='s'
%     image(sign(real(Gamma2Plot)).*sqrt(abs(real(Gamma2Plot)/GammaRealMax))*128+128);
%     ticks=sign(tx).*sqrt(mtix*abs(tx));
%     GammaRealMax=sqrt(GammaRealMax);
% else
%     image(sign(real(Gamma2Plot)).*(abs(real(Gamma2Plot)/GammaRealMax))*128+128);
% end;
% colorbar_rwb(figure1,GammaRealMax,ticks,labels);
% label_boxes(nOrbitals/2,numl,tickx);
% if isunix
%     k = strfind(flnm, '/');
%     if isempty(k)
%             stringp=['/tmp/',flnm(1:length(flnm)-4),'.pdf']
%     else
%         stringp=['/tmp/',flnm(k(length(k))+1:length(flnm)-4),'.pdf']
%     end;
%     %stringp=['/tmp/Gapfunction_realspace_imag.pdf'];
%     print_pdf(stringp);
% end;
% 
% 
% 
%         
