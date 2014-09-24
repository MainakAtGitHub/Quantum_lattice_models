function r=ldos_map(inputfile,ldosflnm)

% inputs
%nOrbitals = 10;
%M = 5;
%ita = .001;
% default for FeSe
sublattice=1;
read_input_file=inputfile;
read_input;
read_input_file
if (~exist('mex','var'))
    mex=true
end;
if (~exist('zGridRange','var'))
    zGridRange = [0 21 4 22]
end;
if (~exist('wannier_filename','var'))
    wannier_filename='wannier_FeSe_4d_matrix_v2.mat'
end;
if (~exist('Greensenergy','var'))
    Greensenergy=0.0084
end;
if (~exist('diagonal_GF','var'))
    diagonal_GF=false;
end;
if ~(exist('singular_quad','var'))
    singular_quad=true;
end;
if ~singular_quad
        sqstring='sum';
end;
E = Greensenergy;
% some double code with impurity_dos (please check, if making
% modifications)
if nargin < 2
    set_ldosfilename;
%LDOSfileName0 = [casestring,'_Vimp_', num2str(Vimp),  '_N_', num2str(N)];
%LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_ita_', num2str(ita),'_e_',num2str(E)];
else
    % input of another filename also accepted
    LDOSfileName=ldosflnm
end;
    
load(LDOSfileName,'-mat')
if diagonal_GF
    diagonal_string='diag';
    latticeGreens=diag(diag(latticeGreens));
else
    diagonal_string='';
end;
%load('./calc/U_0955/LDOS_FeSe_Tom__Vimp_5_N_15_M_10_ita_0.001_e_-0.0084','-mat');
%lattice_greens_supercell_FeSe_N_15_M_9_U_0955_Vimp_5_E_minPt0084.mat
shift = [51 51 41];
%sizeWannier = [101 101 81];
RDiscrete = [40 40 80];
load(wannier_filename,'-mat');
szw=size(wannierValues);
if length(szw)==2
    % Fix for 2D maps
    szw(3)=1;
end;
sizeWannier=szw(1:3);

nBands = size(latticeGreens,1);
if (~exist('N','var'))
    N1 = sqrt(nBands/nOrbitals);
    if ~(N1==N)
        disp(['Error: wrong settings for N:',num2str(N),' vs. ', num2str(N1)]);
    end
end;
if (~exist('xrange','var'))
    % change default behavior: calculate the whole grid
    xrange=N/2
end;
% Local greens function
% Wannier vector
% be a little more general to allow for maps that are not centered at the
% central site
if numel(xrange)==1
    xrange(2)=xrange;
    xrange(1)=-xrange(1);
end;
if (~exist('yrange','var'))
    yrange=xrange;
end;
%xGridRange = -80:60;
if sublattice==1
    xGridRange_limits = [round(RDiscrete(1)/2*(xrange(1)-1)),round(RDiscrete(1)/2*xrange(2))];
    %yGridRange = -60:80;
    yGridRange_limits = [round(RDiscrete(2)/2*yrange(1)),round(RDiscrete(2)/2*(yrange(2)+1))];
elseif sublattice==0
    xGridRange_limits = [round(RDiscrete(1)*xrange(1)),round(RDiscrete(1)*xrange(2))];
    %yGridRange = -60:80;
    yGridRange_limits = [round(RDiscrete(2)*yrange(1)),round(RDiscrete(2)*yrange(2))];
end;
wannierValuesreshape=reshape(wannierValues,sizeWannier(1)*sizeWannier(2)*sizeWannier(3),nOrbitals);
wannierValuesreshape=permute(wannierValuesreshape,[2,1]);
%fileName = ['supercell_local_ldos_FeSe_U_0955_Vimp_5','_N_',num2str(N),'_M_',num2str(M),'_E_',num2str(E),'_ita_',num2str(ita)];
for zGridPoint = zGridRange
    disp(['Calculating ',my_int2str(zGridPoint), 'of (',my_int2str(zGridRange(1)),'..',my_int2str(zGridRange(numel(zGridRange))),')']);
    % calling external function for calculation
    tic
    if ~mex
        localLdos=ldos_map_calc(N,nOrbitals,xGridRange_limits,yGridRange_limits,zGridPoint,RDiscrete,shift,sizeWannier,wannierValuesreshape,latticeGreens);
    else
        localLdos=ldos_map_calc2(N,nOrbitals,xGridRange_limits,yGridRange_limits,zGridPoint,RDiscrete,shift,sizeWannier,wannierValuesreshape,latticeGreens);
        %localLdos=ldos_map_calc_mex(int32(N),int32(nOrbitals),int32(xGridRange_limits),int32(yGridRange_limits),int32(zGridPoint),int32(RDiscrete),int32(shift),int32(sizeWannier),wannierValuesreshape,latticeGreens);
    end;
    toc
    % output of result
    ldosmapfilename=[LDOSfileName,'_z_',num2str(zGridPoint),diagonal_string];
    % save the output together with geometry information necessary to plot
    xGridRange=xGridRange_limits(1):xGridRange_limits(2);
    yGridRange=yGridRange_limits(1):yGridRange_limits(2);
    save(ldosmapfilename,'localLdos','xGridRange','yGridRange','shift','sizeWannier','RDiscrete','sublattice');
end
