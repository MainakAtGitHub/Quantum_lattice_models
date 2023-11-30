function [kxmin,kymin]=plot_banana_1band(inputfile,energy)
% calculate the positions of the ends of the banana shaped equienergylines
% for a 1-band SC (assuming that the result of the homogeneous calculation
% has been written out)

% load the input file to set the variables, gave up the old .mat file
% format
read_input_file=inputfile;
read_input;
read_input_file
M=N*M;
infile=[BdGfileName,'_M_', num2str(M)];
%save([outfile1,'_normal'], 'kSpaceEigenValuesNormal', 'kx', 'ky', '-mat');
load([infile,'_SC'],'-mat');
% cut out first quadrant to speed up calculation and avoid multiples of pi
% in the result
%kx=kx(1:numel(kx)/2);
%ky=ky(1:numel(ky)/2);
% shift by pi
[kxg,kyg]=meshgrid(kx-pi,ky-pi);
%kSpaceEigenValues=kSpaceEigenValues(1:numel(kx),1:numel(kx),:);
kSpaceEigenValues = fftshift(kSpaceEigenValues(:,:,2));
energy=[0.001:0.001:0.033];

% Create figure
figure1 = figure('Position',[200, 50, 200, 150],'PaperUnits','centimeter','PaperPosition',[4 1 12 9]);

% Create axes
axes1 = axes('Parent',figure1);
%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[-0.999285714285714 0.999285714285714]);
%% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-0.999285714285714 0.999285714285714]);
view(axes1,[-25.5 62]);
grid(axes1,'on');

contour3(kxg/pi,kyg/pi,kSpaceEigenValues,energy)
% Create xlabel
xlabel({'k_x/\pi'});

% Create ylabel
ylabel({'k_y/\pi'});

% Create zlabel
zlabel({'E_k'});

% Create colorbar
colorbar('peer',axes1);
[pathstr,name,ext] = fileparts(inputfile) 
print_pdf(['/tmp/banana_',name,'2.pdf'])

view(axes1,[-0.5 90]);
print_pdf(['/tmp/banana_',name,'_top2.pdf'])
