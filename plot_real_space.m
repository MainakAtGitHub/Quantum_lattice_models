function f=plot_real_space(inputfile,scale,tickx,nOrbitals,sublattice)
if nargin < 1
    inputfile='BdG_homogeneous_FeSe_Toms_BS_6Dec13_GammaCut_3_N_9.mat'
end;
if nargin <2
    scale = 's';
end;
if nargin < 3 
    tickx={'$d_{z^2}$','$d_{x^2-y^2}$','$d_{yz}$','$d_{xz}$','$d_{xy}$'}; % orbitals order for FeSe (Tom)
end;
if nargin <4
    nOrbitals=5
end;
if nargin <5
    sublattice=0
end;
fsz=20;
load(inputfile,'-mat'); % BdG_homogeneous_FeSe_Toms_BS_6Dec13_GammaCut_3_N_9.mat
% no need to set the ticks any more
% ticks = -10:5:10;

% extract delta_i0 from delta_ij
N = sqrt(size(delta,1)/nOrbitals);
centerCell = [ceil(N/2) ceil(N/2)];
deltaCenter = zeros(nOrbitals, nOrbitals, N^2);
latticeVectorsDelta = zeros(N^2,2);
count = 0;
for ix = -(ceil(N/2)-1):(ceil(N/2)-1)
    for iy = -(ceil(N/2)-1):(ceil(N/2)-1)
        count = count + 1;
        iCell = centerCell + [ix iy];
        latticeVectorsDelta(count,:) = [ix iy];
        [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, iCell, centerCell);
        deltaCenter(:,:,count) = delta(iRange, jRange);
    end
end

% some output to quantify the homogeneous gap
impCell = [ceil(N/2) ceil(N/2)];
impNNCell = impCell + [0 1];
[iImpNNRange, jImpNNRange] = find_lattice_translation_index(N, nOrbitals, impNNCell, impCell);
if abs(sublattice)==1
    iNNsiteRange = iImpNNRange(1:nOrbitals/2);
    jNNsiteRange = jImpNNRange(1:nOrbitals/2);
    iNNNsiteRange = iImpNNRange((1+nOrbitals/2):nOrbitals);
    jNNNsiteRange = jNNsiteRange;
else
    iNNsiteRange = iImpNNRange(1:nOrbitals);
    jNNsiteRange = jImpNNRange(1:nOrbitals);
    impNNNCell = impCell + [1 1];
    [iImpNNNRange, jImpNNNRange] = find_lattice_translation_index(N, nOrbitals, impNNNCell, impCell);
    iNNNsiteRange = iImpNNNRange(1:nOrbitals);
    jNNNsiteRange = jImpNNNRange(1:nOrbitals);
end;
[iImpNN, jImpNN] = find_lattice_translation_index(N, nOrbitals, impCell, impCell);
deltaOnsite = max(max(abs(delta(iImpNN, jImpNN))));
deltaMaxNN = max(max(abs(delta(iNNsiteRange, jNNsiteRange))));
deltaMaxNNN = max(max(abs(delta(iNNNsiteRange, jNNNsiteRange))));
% give back three numbers that classify the state
f=[deltaOnsite, deltaMaxNN, deltaMaxNNN];
if ~isempty(scale)
delta = deltaCenter;
if ~sublattice==0
%convert 2Fe to 1Fe cell
latticeVectors1Fe = [];
delta1Fe = [];
count = 0;
for i = -(ceil(N/2)-1):(ceil(N/2)-1)
    for j = -(ceil(N/2)-1):(ceil(N/2)-1)
        count = count + 1;
        delta2Fe = delta(:,:,count);
        switch sublattice
            case 1
                % FeSe
                latticeVectors1Fe = [latticeVectors1Fe; i+j, j-i; i+j, j-i+1];
                delta1Fe = [delta1Fe delta2Fe(1:nOrbitals/2,:)];
            case -1
                % LiFeAs (not fixed yet)
             %   latticeVectors1Fe = [latticeVectors1Fe; i+j j-i; i+j j-i-1];
                     %           latticeVectors1Fe = [latticeVectors1Fe; i+j j-i; i+j j-i-1];
                                latticeVectors1Fe = [latticeVectors1Fe; i+j j-i; i+j j-i+1];

                % tmp=delta2Fe((1:nOrbitals/2)+nOrbitals/2,:);
                tmp= delta2Fe(1:nOrbitals/2,:);
                %tmp(:,(1:nOrbitals/2)+nOrbitals/2)=0;
                %tmp(:,(1:nOrbitals/2))=0;
                delta1Fe = [delta1Fe tmp];
        end;
     %   delta1Fe = [delta1Fe delta2Fe(1:nOrbitals/2,:)];
    end
end

% reshape delta
% zoom half way in
plotN=ceil((round(N/2)-1)/2)*2+1;
delta1Fe =  reshape(delta1Fe,nOrbitals/2, nOrbitals/2, 2*N^2);
delta2Plot = zeros(nOrbitals/2*plotN, nOrbitals/2*plotN);
effOrbitals=nOrbitals/2;
else
    plotN=3;
    ceil((round(N/2)-1)/2)*2+1;
    delta1Fe=reshape(delta,nOrbitals,nOrbitals,N^2);
    delta2Plot = zeros(nOrbitals*plotN, nOrbitals*plotN);
    effOrbitals=nOrbitals;
    latticeVectors1Fe=latticeVectorsDelta;
end
for iOrbital = 1:effOrbitals
    for jOrbital = 1:effOrbitals
        deltaBlock = zeros(plotN,plotN);
        for i = 1:plotN
            for j = 1:plotN
                latticeVector = [j-ceil(plotN/2), ceil(plotN/2)-i];
                ind = find((latticeVectors1Fe(:,1) == latticeVector(1)) & (latticeVectors1Fe(:,2) == latticeVector(2)));
                deltaBlock(i,j) = delta1Fe(iOrbital, jOrbital, ind);
            end
        end
        delta2Plot(((iOrbital-1)*plotN + 1):iOrbital*plotN, ((jOrbital-1)*plotN + 1):jOrbital*plotN ) = deltaBlock;
    end
end
% convert to meV
delta2Plot = 1000*delta2Plot;
%r=realspaceplot(data2plot,numl,tickx,flnm,scale)
cptn='\Delta_{RR''}^{\mu\nu} [meV]';
[~,h,stringp]=realspaceplot(delta2Plot,plotN,tickx,inputfile,scale,cptn);
if isunix
    % save delta in format as Gamma
    nup1=nUp(1:nOrbitals);
    ndown1=nDown(1:nOrbitals);
    strng=[stringp(1:length(stringp)-4),'_delta_hom.mat',]
    save(strng,'delta1Fe','latticeVectors1Fe','nup1','ndown1','mu');
end;
end
% figure1=figure
% 
% global colorred
% colorred=false; % set colorscale
% global fsz
% fsz = 25; % font
%  deltaRealMax = max(abs(real(delta2Plot(:)))); % for setting colormap
% % if scale == 's'
% %     image(sign(real(delta2Plot)).*sqrt(abs(real(delta2Plot)/deltaRealMax))*128+128);
% % else
% %     image(sign(real(delta2Plot)).*(abs(real(delta2Plot)/deltaRealMax))*128+128);
% % end
% % colorbar_rwb(figure1,deltaRealMax,ticks);
% % label_boxes(nOrbitals/2,numl,tickx);
% 
% % set the scale
% lm=log(deltaRealMax)/log(10);
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
%     image(sign(real(delta2Plot)).*sqrt(abs(real(delta2Plot)/deltaRealMax))*128+128);
%     ticks=sign(tx).*sqrt(mtix*abs(tx));
%     deltaRealMax=sqrt(deltaRealMax);
% else
%     image(sign(real(delta2Plot)).*(abs(real(delta2Plot)/deltaRealMax))*128+128);
% end;
% colorbar_rwb(figure1,deltaRealMax,ticks,labels);
% label_boxes(nOrbitals/2,numl,tickx);
% if isunix
%     k = strfind(inputfile, '/');
%     if isempty(k)
%             stringp=['/tmp/',inputfile(1:length(inputfile)-4),'.pdf']
%     else
%         stringp=['/tmp/',inputfile(k(length(k))+1:length(inputfile)-4),'.pdf']
%     end;
%     %stringp=['/tmp/Gapfunction_realspace_imag.pdf'];
%     print_pdf(stringp);
% end;