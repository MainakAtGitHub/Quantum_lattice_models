function f=plot_real_space(inputfile,sqrtscale,tickx,nOrbitals)
if nargin < 1
    inputfile='BdG_homogeneous_FeSe_Toms_BS_6Dec13_GammaCut_3_N_9.mat'
end;
if nargin <2
    sqrtscale=true;
end;
if nargin < 3 
    tickx={'$d_{z^2}$','$d_{x^2-y^2}$','$d_{yz}$','$d_{xz}$','$d_{xy}$'}; % orbitals order
end;
if nargin <4
    nOrbitals=10;
end;
load(inputfile); % BdG_homogeneous_FeSe_Toms_BS_6Dec13_GammaCut_3_N_9.mat
ticks = -10:5:10;

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
        latticeVectors1Fe = [latticeVectors1Fe; i+j j-i; i+j j-i+1];
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
global colorred
colorred=false; % set colorscale
figure1 = figure;
global fsz
fsz = 25; % font
deltaRealMax = max(abs(real(delta2Plot(:)))); % for setting colormap
scale = 's';
if scale == 's'
    image(sign(real(delta2Plot)).*sqrt(abs(real(delta2Plot)/deltaRealMax))*128+128);
else
    image(sign(real(delta2Plot)).*(abs(real(delta2Plot)/deltaRealMax))*128+128);
end
colorbar_rwb(figure1,deltaRealMax,ticks);
label_boxes(5,numl,tickx);

