function [max_gap,min_gap]=plot_delta_map(inputfile,which_del)

% Modified homogeneous_dos.m
% takes \Delta_ij as input and constructs \Delta_i0.
% which_del can take values 1,2,3,4 corresponding to DeltaUpUp,
% DeltaUpDown, DeltaDownUp, DeltaDownDown respectively

if nargin <2 
    which_del=0;
end
if nargin <1
    % load relevant files
    TB_file='TB_hamiltonian_FeSe_2D.mat'
    BdGfileName='BdG_homogeneous_FeSe_Milan_GammaCut_2_N_9(1).mat'
    Gamma_file='Gamma_FeSe_Milan_GammaCut_2.mat'
    M = input('Enter no of k-points   ');% no of K points in x
    ita = input('Enter ita   '); % broadening
    firstEnergy = input('Enter starting energy   ');
    lastEnergy = input('Enter last energy   ');
    nEnergyPoints = input('Enter no of energy points   ');
else
    % load the input file to set the variables, gave up the old .mat file
    % format
     read_input_file=inputfile;
     read_input;
     read_input_file
end;
if (~exist('Vimp','var'))
    if abs(Vimp)>0
        disp('Warning: finite impurity potential, not homogeneous case.')
    end;
end;
if ~(exist('singular_quad','var'))
    singular_quad=true;
end;
if ~(exist('write_states','var'))
    write_states=false;
end;
if ~(exist('calc_dos','var'))
    calc_dos=true;
end;
if ~singular_quad
        sqstring='sum';
else
    sqstring='';
end;
load(TB_file,'-mat');

if ~exist('nOrbitals','var') % if {~exist('nOrbitals','var')...end} added on Jan2021......because nOrbitals already exists in inputfile as a keyword in newer runs
nOrbitals = size(TBparameters,1);
end

load(Gamma_file,'-mat');
load(BdGfileName,'-mat');


if which_del==1
    delta=delta(1:size(delta,1)/2, 1:size(delta,1)/2);
elseif which_del==2
    delta=delta(1:size(delta,1)/2, (size(delta,1)/2 + 1):(size(delta,1)/2 + size(delta,1)/2));
elseif which_del==3
    delta=(delta((size(delta,1)/2 + 1):(size(delta,1)/2 + size(delta,1)/2), 1:size(delta,1)/2));
elseif which_del==4
    delta=delta((size(delta,1)/2 + 1):(size(delta,1)/2 + size(delta,1)/2), (size(delta,1)/2 + 1):(size(delta,1)/2 + size(delta,1)/2));
end


N = sqrt(size(delta,1)/nOrbitals);
if exist('latticeVectorsSC','var')
    nUnitCellsDelta = size(latticeVectorsSC,1);
else
    latticeVectorsSC=Gammafull.latt;
    nUnitCellsDelta = size(Gammafull.latt,1);
end


deltamap = zeros(N, N, nOrbitals,nOrbitals);
% read out delta and store it into map of matrices
% to do generalize to non-rectangular systems
for N1=1:N
    for N2=1:N
        jCell = [N1 N2];
        for i = 1:nUnitCellsDelta
            iCell = jCell + latticeVectorsSC(i,:);
            cellvector=periodic_latticevectors(iCell,N);
            [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, cellvector, jCell);
            position=jCell+latticeVectorsSC(i,:)/2;
            % now add the contribution to the 4 neighbored lattice points
            pos(1,:)=periodic_latticevectors(floor(position),N);
            pos(2,:)=periodic_latticevectors(ceil(position),N);
            pos(3,:)=periodic_latticevectors([ceil(position(1)),floor(position(2))],N);
            pos(4,:)=periodic_latticevectors([floor(position(1)),ceil(position(2))],N);
            for ps=1:4
                deltamap(pos(ps,1),pos(ps,2),:,:) =deltamap(pos(ps,1),pos(ps,2),:,:)+shiftdim(delta(iRange, jRange).*conj(delta(iRange, jRange)),-2);
            end
        end
    end
end

% Mainak (adding horizontal and vertical separate maps)##############
deltamapH = zeros(N, N, nOrbitals,nOrbitals);
% read out delta and store it into map of matrices
% to do generalize to non-rectangular systems
for N1=1:N
    for N2=1:N
        jCell = [N1 N2];
        for i = [1,5]
            iCell = jCell + latticeVectorsSC(i,:);
            cellvector=periodic_latticevectors(iCell,N);
            [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, cellvector, jCell);
            position=jCell+latticeVectorsSC(i,:)/2;
            % now add the contribution to the 4 neighbored lattice points
            pos(1,:)=periodic_latticevectors(floor(position),N);
            pos(2,:)=periodic_latticevectors(ceil(position),N);
            pos(3,:)=periodic_latticevectors([ceil(position(1)),floor(position(2))],N);
            pos(4,:)=periodic_latticevectors([floor(position(1)),ceil(position(2))],N);
            for ps=1:4
                deltamapH(pos(ps,1),pos(ps,2),:,:) =deltamapH(pos(ps,1),pos(ps,2),:,:)+shiftdim(delta(iRange, jRange),-2);
            end
        end
    end
end
% Mainak ####################################################

% Mainak (adding horizontal and vertical separate maps)##############
deltamapV = zeros(N, N, nOrbitals,nOrbitals);
% read out delta and store it into map of matrices
% to do generalize to non-rectangular systems
for N1=1:N
    for N2=1:N
        jCell = [N1 N2];
        for i = [2,4]
            iCell = jCell + latticeVectorsSC(i,:);
            cellvector=periodic_latticevectors(iCell,N);
            [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, cellvector, jCell);
            position=jCell+latticeVectorsSC(i,:)/2;
            % now add the contribution to the 4 neighbored lattice points
            pos(1,:)=periodic_latticevectors(floor(position),N);
            pos(2,:)=periodic_latticevectors(ceil(position),N);
            pos(3,:)=periodic_latticevectors([ceil(position(1)),floor(position(2))],N);
            pos(4,:)=periodic_latticevectors([floor(position(1)),ceil(position(2))],N);
            for ps=1:4
                deltamapV(pos(ps,1),pos(ps,2),:,:) =deltamapV(pos(ps,1),pos(ps,2),:,:)+shiftdim(delta(iRange, jRange),-2);
            end
        end
    end
end
% Mainak ####################################################

delta_map_single=zeros(N, N);
% do some plotting of the maps
for N1=1:N
    for N2=1:N
        delta_map_single(N1,N2)=sum(sum(squeeze(deltamap(N1,N2,:,:))));
    end
end
delta_map_single=sqrt(delta_map_single);
figure
imagesc(delta_map_single');
% give back maximum and minimum for purpose of plotting phase diagrams
max_gap=max(delta_map_single(:));
min_gap=min(delta_map_single(:));
axis square;
colorbar
[filepath,name,ext]=fileparts(BdGfileName);
% also load the impurity positions and plot
if ~isnumeric(Vimp)
load(Vimp);
hold on
impCell = [ceil(N/2) ceil(N/2)];
for s=1:size(imp_vec,1)
    imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
end
plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
end
print_pdf([filepath,filesep,name,'_gapmap.pdf']);

% Mainak (plotting hor gaps)#######################
delta_map_singleH=zeros(N, N);
% do some plotting of the maps
for N1=1:N
    for N2=1:N
        delta_map_singleH(N1,N2)=sum(sum(squeeze(deltamapH(N1,N2,:,:))));
    end
end
%%%delta_map_singleH=sqrt(delta_map_singleH);
figure
imagesc(delta_map_singleH');
axis square;
colorbar
[filepath,name,ext]=fileparts(BdGfileName);
% also load the impurity positions and plot
if ~isnumeric(Vimp)
load(Vimp);
hold on
impCell = [ceil(N/2) ceil(N/2)];
for s=1:size(imp_vec,1)
    imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
end
plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
end
print_pdf([filepath,filesep,name,'_gapmap_hor.pdf']);
% Mainak###################################################

% Mainak (plotting ver gaps)#######################
delta_map_singleV=zeros(N, N);
% do some plotting of the maps
for N1=1:N
    for N2=1:N
        delta_map_singleV(N1,N2)=sum(sum(squeeze(deltamapV(N1,N2,:,:))));
    end
end
%%%delta_map_singleV=sqrt(delta_map_singleV);
figure
imagesc(delta_map_singleV');
axis square;
colorbar
[filepath,name,ext]=fileparts(BdGfileName);
% also load the impurity positions and plot
if ~isnumeric(Vimp)
load(Vimp);
hold on
impCell = [ceil(N/2) ceil(N/2)];
for s=1:size(imp_vec,1)
    imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
end
plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
end
print_pdf([filepath,filesep,name,'_gapmap_ver.pdf']);
% Mainak###################################################