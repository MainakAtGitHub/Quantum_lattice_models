function [max_gap,min_gap]=plot_delta_map(inputfile,which_del,deltachannel,specific_orb,flash_fig)
if ~exist('flash_fig','var')
    flash_fig='on'; %%%%%% flash_fig = 'on' or 'off'
end
% global save_d_ord_here1; % temporary global var for dirty dwave order par 
% global noise_check_dwave3;
% global bond_ord_hor_ver;
% Modified homogeneous_dos.m
% takes \Delta_ij as input and constructs \Delta_i0.
% which_del can take values 1,2,3,4 corresponding to DeltaUpUp,
% DeltaUpDown, DeltaDownUp, DeltaDownDown respectively

if nargin < 4
    specific_orb=0;
end

if nargin <3 
    deltachannel=0; % setting default to be the singlet channel
end
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

if ~exist('NNN_pairing','var')
    NNN_only_pairing=false;
end

if ~exist('spin_and_nambu','var')
    spin_and_nambu=false;
end
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
if exist('pos_file','var')
    load(pos_file,'-mat');
else
    if numel(N)==1
        N = [N,N];
    end
    [tempxx,tempyy]=meshgrid(-N(1)/2+0.5*mod(N(1)+1,2):N(1)/2-0.5*mod(N(1)+1,2),-N(2)/2+0.5*mod(N(2)+1,2):N(2)/2-0.5*mod(N(2)+1,2));
    r=[tempxx(:),tempyy(:)];
    clear tempxx tempyy;
end
if spin_and_nambu
    if which_del==1
        delta=delta(1:size(delta,1)/2, 1:size(delta,1)/2);
    elseif which_del==2
        delta=delta(1:size(delta,1)/2, (size(delta,1)/2 + 1):(size(delta,1)/2 + size(delta,1)/2));
    elseif which_del==3
        delta=(delta((size(delta,1)/2 + 1):(size(delta,1)/2 + size(delta,1)/2), 1:size(delta,1)/2));
    elseif which_del==4
        delta=delta((size(delta,1)/2 + 1):(size(delta,1)/2 + size(delta,1)/2), (size(delta,1)/2 + 1):(size(delta,1)/2 + size(delta,1)/2));
    end
    
    if any(deltachannel==[1,2,3])
        if exist('tripletOrderFile','var')
            load(tripletOrderFile);
            tripletOrderMatrix=lattice_translation(N, tripletOrder, tripletOrderVecs);
            delta=tripletOrderMatrix.*delta;
        end
    end
end
% N = sqrt(size(delta,1)/nOrbitals);
if numel(N)==1
    N=[N,N];
end
if exist('latticeVectorsSC','var')
    nUnitCellsDelta = size(latticeVectorsSC,1);
else
    latticeVectorsSC=Gammafull.latt;
    nUnitCellsDelta = size(Gammafull.latt,1);
end


deltamap = zeros(N(1), N(2), nOrbitals,nOrbitals);
% read out delta and store it into map of matrices
% to do generalize to non-rectangular systems
for N1=1:N(1)
    for N2=1:N(2)
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
%               deltamap(pos(ps,1),pos(ps,2),:,:) =deltamap(pos(ps,1),pos(ps,2),:,:)+shiftdim(delta(iRange, jRange).*conj(delta(iRange, jRange)),-2);
            end
        end
    end
end

% Mainak (adding horizontal and vertical separate maps)##############

if nUnitCellsDelta > 5
    
    
    deltamapH = zeros(N(1), N(2), nOrbitals,nOrbitals);
    % read out delta and store it into map of matrices
    % to do generalize to non-rectangular systems
    for N1=1:N(1)
        for N2=1:N(2)
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
    deltamapH=deltamapH/8; % Normalization
    
    deltamapV = zeros(N(1), N(2), nOrbitals,nOrbitals);
    % read out delta and store it into map of matrices
    % to do generalize to non-rectangular systems
    for N1=1:N(1)
        for N2=1:N(2)
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
    deltamapV=deltamapV/8;
    
    deltamapD_xpy = zeros(N(1), N(2), nOrbitals,nOrbitals);
    % read out delta and store it into map of matrices
    % to do generalize to non-rectangular systems
    for N1=1:N(1)
        for N2=1:N(2)
            jCell = [N1 N2];
            for i = [6,9]
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
                    deltamapD_xpy(pos(ps,1),pos(ps,2),:,:) =deltamapD_xpy(pos(ps,1),pos(ps,2),:,:)+shiftdim(delta(iRange, jRange),-2);
                end
            end
        end
    end
    % Mainak ####################################################
    deltamapD_xpy=deltamapD_xpy/8;
    
    deltamapD_xmy = zeros(N(1), N(2), nOrbitals,nOrbitals);
    % read out delta and store it into map of matrices
    % to do generalize to non-rectangular systems
    for N1=1:N(1)
        for N2=1:N(2)
            jCell = [N1 N2];
            for i = [7,8]
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
                    deltamapD_xmy(pos(ps,1),pos(ps,2),:,:) =deltamapD_xmy(pos(ps,1),pos(ps,2),:,:)+shiftdim(delta(iRange, jRange),-2);
                end
            end
        end
    end
    % Mainak ####################################################
    deltamapD_xmy=deltamapD_xmy/8;    
    
    
end

if nUnitCellsDelta <= 5 %%%% restricting to the NN only or NNN only pairing
    
    deltamapH = zeros(N(1), N(2), nOrbitals,nOrbitals);
    % read out delta and store it into map of matrices
    % to do generalize to non-rectangular systems
    for N1=1:N(1)
        for N2=1:N(2)
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
    deltamapH=deltamapH/8; % Normalization
end

% Mainak (adding horizontal and vertical separate maps)##############

if nUnitCellsDelta <= 5
    deltamapV = zeros(N(1), N(2), nOrbitals,nOrbitals);
    % read out delta and store it into map of matrices
    % to do generalize to non-rectangular systems
    for N1=1:N(1)
        for N2=1:N(2)
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
    deltamapV=deltamapV/8;
end

delta_map_single=zeros(N(1), N(2));
% do some plotting of the maps
for N1=1:N(1)
    for N2=1:N(2)
        if specific_orb
            delta_map_single(N1,N2)=(((deltamap(N1,N2,specific_orb,specific_orb))));
        else
            delta_map_single(N1,N2)=sum(sum(squeeze(deltamap(N1,N2,:,:))));
        end
    end
end
delta_map_single=sqrt(delta_map_single);
h=figure; 
set(h, 'Visible', flash_fig); %%%%%%%% uncomment if want to see the figure
pcolor(delta_map_single');
title(['rms NN gap',', system avg. rms NN gap = ',num2str(sum(sum(delta_map_single'))/numel(delta_map_single'))]);
% give back maximum and minimum for purpose of plotting phase diagrams
max_gap=max(delta_map_single(:));
min_gap=min(delta_map_single(:));
axis equal;
colorbar
[filepath,name,ext]=fileparts(BdGfileName);
% also load the impurity positions and plot
if ~isnumeric(Vimp)
load(Vimp);
hold on
impCell = [ceil(N(1)/2) ceil(N(2)/2)];
for s=1:size(imp_vec,1)
    imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
end
plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
end
print_pdf([filepath,filesep,name,'_gapmap.pdf']);



if nUnitCellsDelta > 5 %%%%%%%%%%%%%%%%%%%% plotting NN and NNN
   
    % Mainak (plotting hor gaps)#######################
    delta_map_singleH=zeros(N(1), N(2));
    % do some plotting of the maps
    for N1=1:N(1)
        for N2=1:N(2)
            if specific_orb
                delta_map_singleH(N1,N2)=(((deltamapH(N1,N2,specific_orb,specific_orb))));
            else
                delta_map_singleH(N1,N2)=sum(sum(squeeze(deltamapH(N1,N2,:,:))));
            end
        end
    end
    %%%delta_map_singleH=sqrt(delta_map_singleH);
    h=figure;
    set(h, 'Visible', flash_fig);
    pcolor(abs(delta_map_singleH'));
    axis equal;
    colorbar;
    hold on;
    dummyreal=real(delta_map_singleH');
    dummyimag=imag(delta_map_singleH');
    quiver(dummyreal,dummyimag);
%     if NNN_only_pairing==true
%         title('NNN(x+y) SC gap');
%     else
    title('Horizontal SC gap');
%     end
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end        
    [filepath,name,ext]=fileparts(BdGfileName);
    print_pdf([filepath,filesep,name,'_gapmap_hor.pdf']);
    
    
    h=figure;
    set(h, 'Visible', flash_fig);
    subplot(1,2,1);
    pcolor(real(delta_map_singleH'));
    axis equal;
    colorbar
    subplot(1,2,2);
    pcolor(imag(delta_map_singleH'));
    axis equal;
    colorbar
    [filepath,name,ext]=fileparts(BdGfileName);
    % also load the impurity positions and plot
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end    
    
    %%%%%%%%%%%%%%5 ver gap below
    delta_map_singleV=zeros(N(1), N(2));
    % do some plotting of the maps
    for N1=1:N(1)
        for N2=1:N(2)
            if specific_orb
                delta_map_singleV(N1,N2)=(((deltamapV(N1,N2,specific_orb,specific_orb))));
            else
                delta_map_singleV(N1,N2)=sum(sum(squeeze(deltamapV(N1,N2,:,:))));
            end
        end
    end
    %%%delta_map_singleV=sqrt(delta_map_singleV);
    h=figure;
    set(h, 'Visible', flash_fig);
    pcolor(abs(delta_map_singleV'));
    axis equal;
    colorbar;
    hold on;
    dummyreal=real(delta_map_singleV');
    dummyimag=imag(delta_map_singleV');
    quiver(dummyreal,dummyimag);
%     if NNN_only_pairing==true
%         title('NNN(x-y) SC gap');
%     else
    title('Vertical SC gap');
%     end
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end        
    [filepath,name,ext]=fileparts(BdGfileName);
    print_pdf([filepath,filesep,name,'_gapmap_ver.pdf']);
    
    h=figure;
    set(h, 'Visible', flash_fig);
    subplot(1,2,1);
    pcolor(real(delta_map_singleV'));
    axis equal;
    colorbar
    subplot(1,2,2);
    pcolor(imag(delta_map_singleV'));
    axis equal;
    colorbar
    
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end     
    
    %%%%%%%%%%%%%%%%%%%%diagonal xpy plotting below
    % #######################
    delta_map_singleD_xpy=zeros(N(1), N(2));
    % do some plotting of the maps
    for N1=1:N(1)
        for N2=1:N(2)
            if specific_orb
                delta_map_singleD_xpy(N1,N2)=(((deltamapD_xpy(N1,N2,specific_orb,specific_orb))));
            else
                delta_map_singleD_xpy(N1,N2)=sum(sum(squeeze(deltamapD_xpy(N1,N2,:,:))));
            end
        end
    end
    %%%delta_map_singleH=sqrt(delta_map_singleH);
    h=figure;
    set(h, 'Visible', flash_fig);
    dmdlt=delta_map_singleD_xpy';
    [dltpltx,dltplty]=meshgrid(1:N(1),1:N(2));
    scatter(dltpltx(:),dltplty(:),150,abs(dmdlt(:)),'filled');
    axis equal;
    colorbar;
    hold on;
    dummyreal=real(delta_map_singleD_xpy');
    dummyimag=imag(delta_map_singleD_xpy');
    quiver(dummyreal,dummyimag);
%     if NNN_only_pairing==true
%         title('NNN(x+y) SC gap');
%     else
    title('NNN (x+y) SC gap');
%     end
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end        
    [filepath,name,ext]=fileparts(BdGfileName);
    print_pdf([filepath,filesep,name,'_gapmap_(x+y).pdf']);
    
    
    h=figure;
    set(h, 'Visible', flash_fig);
    subplot(1,2,1);
    pcolor(real(delta_map_singleD_xpy'));
    axis equal;
    colorbar
    subplot(1,2,2);
    pcolor(imag(delta_map_singleD_xpy'));
    axis equal;
    colorbar
    [filepath,name,ext]=fileparts(BdGfileName);
    % also load the impurity positions and plot
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end   
    
    
    %%%%%%%%%%%%%%%%%%%%diagonal xmy plotting below
    % #######################
    delta_map_singleD_xmy=zeros(N(1), N(2));
    % do some plotting of the maps
    for N1=1:N(1)
        for N2=1:N(2)
            if specific_orb
                delta_map_singleD_xmy(N1,N2)=(((deltamapD_xmy(N1,N2,specific_orb,specific_orb))));
            else
                delta_map_singleD_xmy(N1,N2)=sum(sum(squeeze(deltamapD_xmy(N1,N2,:,:))));
            end
        end
    end
    %%%delta_map_singleH=sqrt(delta_map_singleH);
    h=figure;
    set(h, 'Visible', flash_fig);
    pcolor(abs(delta_map_singleD_xmy'));
    axis equal;
    colorbar;
    hold on;
    dummyreal=real(delta_map_singleD_xmy');
    dummyimag=imag(delta_map_singleD_xmy');
    quiver(dummyreal,dummyimag);
%     if NNN_only_pairing==true
%         title('NNN(x+y) SC gap');
%     else
    title('NNN (x-y) SC gap');
%     end
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end        
    [filepath,name,ext]=fileparts(BdGfileName);
    print_pdf([filepath,filesep,name,'_gapmap_(x-y).pdf']);
    
    
    h=figure;
    set(h, 'Visible', flash_fig);
    subplot(1,2,1);
    pcolor(real(delta_map_singleD_xmy'));
    axis equal;
    colorbar
    subplot(1,2,2);
    pcolor(imag(delta_map_singleD_xmy'));
    axis equal;
    colorbar
    [filepath,name,ext]=fileparts(BdGfileName);
    % also load the impurity positions and plot
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end     
    
    
    
end             %%%%%%%%%%%%%%%%%%%% plotting NN and NNN

if nUnitCellsDelta <= 5
    % Mainak (plotting hor gaps)#######################
    delta_map_singleH=zeros(N(1), N(2));
    % do some plotting of the maps
    for N1=1:N(1)
        for N2=1:N(2)
            if specific_orb
                delta_map_singleH(N1,N2)=(((deltamapH(N1,N2,specific_orb,specific_orb))));
            else
                delta_map_singleH(N1,N2)=sum(sum(squeeze(deltamapH(N1,N2,:,:))));
            end
        end
    end
    %%%delta_map_singleH=sqrt(delta_map_singleH);
    h=figure;
    set(h, 'Visible', flash_fig);
    pcolor(abs(delta_map_singleH'));
    axis equal;
    colorbar;
    hold on;
    dummyreal=real(delta_map_singleH');
    dummyimag=imag(delta_map_singleH');
    quiver(dummyreal,dummyimag);
    if NNN_only_pairing==true
        title('NNN(x+y) SC gap');
    else
        title('Horizontal SC gap');
    end

    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end
    
    
    [filepath,name,ext]=fileparts(BdGfileName);
    if NNN_only_pairing==true
        print_pdf([filepath,filesep,name,'_gapmap_(x+y).pdf']);
    else
        print_pdf([filepath,filesep,name,'_gapmap_hor.pdf']);
    end
    
    h=figure;
    set(h, 'Visible', flash_fig);
    subplot(1,2,1);
    pcolor(real(delta_map_singleH'));
    axis equal;
    colorbar
    subplot(1,2,2);
    pcolor(imag(delta_map_singleH'));
    axis equal;
    colorbar
    [filepath,name,ext]=fileparts(BdGfileName);
    % also load the impurity positions and plot
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end
end
% print_pdf([filepath,filesep,name,'_gapmap_hor.pdf']);
% Mainak###################################################

% Mainak (plotting ver gaps)#######################

if nUnitCellsDelta <= 5
    delta_map_singleV=zeros(N(1), N(2));
    % do some plotting of the maps
    for N1=1:N(1)
        for N2=1:N(2)
            if specific_orb
                delta_map_singleV(N1,N2)=(((deltamapV(N1,N2,specific_orb,specific_orb))));
            else
                delta_map_singleV(N1,N2)=sum(sum(squeeze(deltamapV(N1,N2,:,:))));
            end
        end
    end
    %%%delta_map_singleV=sqrt(delta_map_singleV);
    h=figure;
    set(h, 'Visible', flash_fig);
    pcolor(abs(delta_map_singleV'));
    axis equal;
    colorbar;
    hold on;
    dummyreal=real(delta_map_singleV');
    dummyimag=imag(delta_map_singleV');
    quiver(dummyreal,dummyimag);
    if NNN_only_pairing==true
        title('NNN(x-y) SC gap');
    else
        title('Vertical SC gap');
    end
    
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end
        
    [filepath,name,ext]=fileparts(BdGfileName);
    
    if NNN_only_pairing==true
        print_pdf([filepath,filesep,name,'_gapmap_(x-y).pdf']);
    else
        print_pdf([filepath,filesep,name,'_gapmap_ver.pdf']);
    end
    
    h=figure;
    set(h, 'Visible', flash_fig);
    subplot(1,2,1);
    pcolor(real(delta_map_singleV'));
    axis equal;
    colorbar
    subplot(1,2,2);
    pcolor(imag(delta_map_singleV'));
    axis equal;
    colorbar
    
    [filepath,name,ext]=fileparts(BdGfileName);
    
    %%%%%%%%%%%%%%%%%% d-wave ord. par.
    h=figure;
    set(h, 'Visible', flash_fig);
    imagesc((delta_map_singleH'-delta_map_singleV'));
    axis equal;
    colorbar;
    title(['d-wave gap (hor. gap - ver. gap)',' syst. avg. = ',num2str(sum(sum(delta_map_singleH'-delta_map_singleV'))/numel(delta_map_single'))]);
    [filepath,name,ext]=fileparts(BdGfileName);
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end    
    print_pdf([filepath,filesep,name,'_gapmap_d-wave_HorMinusVer.pdf']);
    bond_ord_hor_ver(:,:,1) = delta_map_singleH'; 
    bond_ord_hor_ver(:,:,2) = delta_map_singleV';
    save([filepath,filesep,'bond_ord_hor_ver.mat'],'bond_ord_hor_ver');
%     noise_check_dwave_delta_map=delta_map_singleH'-delta_map_singleV';
%     load(['/home/UFAD/mainak.pal/Desktop/2021stff/Sep2021/more_impurity_configuration_prep_for_hpc/systematic_impurity_calcs_generated/imp_conc_0.15_noise_check_dwave/','noise_check_dwave.mat'],'noise_check_dwave3');
%     noise_check_dwave3 = [noise_check_dwave3;noise_check_dwave_delta_map(:)];
%     save(['/home/UFAD/mainak.pal/Desktop/2021stff/Sep2021/more_impurity_configuration_prep_for_hpc/systematic_impurity_calcs_generated/imp_conc_0.15_noise_check_dwave/','noise_check_dwave.mat'],'noise_check_dwave3');
    % tempsave_d_ord_here1=[save_d_ord_here1;[2*(0.85-n0),sum(sum(delta_map_singleH'-delta_map_singleV'))/numel(delta_map_single')]];
    %

    %%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%% d-wave ord. par. abs
    h=figure;
    set(h, 'Visible', flash_fig);
    imagesc(abs(delta_map_singleH'-delta_map_singleV'));
    axis equal;
    colorbar;
    title(['abs d-wave gap (hor. gap - ver. gap)',' syst. avg. = ',num2str(sum(sum(abs(delta_map_singleH'-delta_map_singleV')))/numel(delta_map_single'))]);
    [filepath,name,ext]=fileparts(BdGfileName);
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end    
    print_pdf([filepath,filesep,name,'_abs_gapmap_d-wave_HorMinusVer.pdf']);
    %

    %%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%% d-wave ord. par. real part
    h=figure;
    set(h, 'Visible', flash_fig);
    imagesc(real(delta_map_singleH'-delta_map_singleV'));
    axis equal;
    colorbar;
    title(['real d-wave gap (hor. gap - ver. gap)',' syst. avg. = ',num2str(sum(sum(real(delta_map_singleH'-delta_map_singleV')))/numel(delta_map_single'))]);
    [filepath,name,ext]=fileparts(BdGfileName);
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end    
    print_pdf([filepath,filesep,name,'_real_gapmap_d-wave_HorMinusVer.pdf']);
    
    %

    %%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%% d-wave ord. par. imag part
    h=figure;
    set(h, 'Visible', flash_fig);
    imagesc(imag(delta_map_singleH'-delta_map_singleV'));
    axis equal;
    colorbar;
    title(['imag d-wave gap (hor. gap - ver. gap)',' syst. avg. = ',num2str(sum(sum(imag(delta_map_singleH'-delta_map_singleV')))/numel(delta_map_single'))]);
    [filepath,name,ext]=fileparts(BdGfileName);
    if ~isnumeric(Vimp)
        load(Vimp);
        hold on
        impCell = [ceil(N(1)/2) ceil(N(2)/2)];
        for s=1:size(imp_vec,1)
            imp_vec_shift(s,:)=periodic_latticevectors(impCell+imp_vec(s,:),N);
        end
        plot(imp_vec_shift(:,1),imp_vec_shift(:,2),'xr');
    end    
    print_pdf([filepath,filesep,name,'_imag_gapmap_d-wave_HorMinusVer.pdf']);
    %%%%%%%%%%%%%%%%%%%
    
    % also load the impurity positions and plot

end
% print_pdf([filepath,filesep,name,'_gapmap_ver.pdf']);
% Mainak###################################################