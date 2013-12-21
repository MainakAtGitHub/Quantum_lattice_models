function f=plot_real_space(inputfile,scale,tickx,nOrbitals)
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
    nOrbitals=10;
end;
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
delta = deltaCenter;

%convert 2Fe to 1Fe cell
latticeVectors1Fe = [];
delta1Fe = [];
count = 0;
for i = -(ceil(N/2)-1):(ceil(N/2)-1)
    for j = -(ceil(N/2)-1):(ceil(N/2)-1)
        count = count + 1;
        delta2Fe = delta(:,:,count);
        % FeSe
         latticeVectors1Fe = [latticeVectors1Fe; i+j j-i; i+j j-i+1];
        % LiFeAs (not fixed yet)
        %latticeVectors1Fe = [latticeVectors1Fe; -j+i i+j; i-j-1 i+j];        
        delta1Fe = [delta1Fe delta2Fe(1:nOrbitals/2,:)];
    end
end

% reshape delta
delta1Fe =  reshape(delta1Fe,nOrbitals/2, nOrbitals/2, 2*N^2);
delta2Plot = zeros(nOrbitals/2*N, nOrbitals/2*N);
for iOrbital = 1:nOrbitals/2
    for jOrbital = 1:nOrbitals/2
        deltaBlock = zeros(N,N);
        for i = 1:N
            for j = 1:N
                latticeVector = [j-ceil(N/2), ceil(N/2)-i];
                ind = find((latticeVectors1Fe(:,1) == latticeVector(1)) & (latticeVectors1Fe(:,2) == latticeVector(2)));
                deltaBlock(i,j) = delta1Fe(iOrbital, jOrbital, ind);
            end
        end
        delta2Plot(((iOrbital-1)*N + 1):iOrbital*N, ((jOrbital-1)*N + 1):jOrbital*N ) = deltaBlock;
    end
end
delta2Plot = 1000*delta2Plot;
numl = N;

%r=realspaceplot(data2plot,numl,tickx,flnm,scale)
r=realspaceplot(delta2Plot,numl,tickx,inputfile,'s');
max(abs(delta2Plot(:)))
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



