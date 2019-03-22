function r=ldos_map(inputfile,ldosflnm)
% derived from ldos_map.m r47
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
if (~exist('nOrbitals','var'))
    nOrbitals = 1
end;
if (~exist('N','var'))
    N = 11
end;
if (~exist('wannier_filename','var'))
    wannier_filename='wannier_FeSe_4d_matrix_v2.mat'
end;
if ~(exist('Greensenergy','var'))
    disp('No Greensenergy given, setting to 0.');
    E=0;
else
    if ~ischar(Greensenergy)
        E=Greensenergy;
    else
        % load the list of energies to be calculated from the given
        % file
        load(Greensenergy,'E');
        le=numel(E);
    end;
end
le=numel(E);
if (~exist('diagonal_GF','var'))
    diagonal_GF=false;
end;
if ~(exist('singular_quad','var'))
    singular_quad=true;
end;
if ~singular_quad
        sqstring='sum';
end;
% E = Greensenergy;
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
    if (diagonal_GF==2)
        diagonal_string='diag_orb'
        % make greens function diagonal in orbitals
        latticeGreens=latticeGreens.*repmat(eye(nOrbitals),N^2,N^2);
    elseif (diagonal_GF==3)
        diagonal_string='diag_local'
        latticeGreens=latticeGreens.*kron(eye(N^2),ones(nOrbitals));
    else
        diagonal_string='diag'
        latticeGreens=diag(diag(latticeGreens));
    end;
else
    diagonal_string='';
end;
%load('./calc/U_0955/LDOS_FeSe_Tom__Vimp_5_N_15_M_10_ita_0.001_e_-0.0084','-mat');
%lattice_greens_supercell_FeSe_N_15_M_9_U_0955_Vimp_5_E_minPt0084.mat
shift = [51 51 41];
%sizeWannier = [101 101 81];
RDiscrete = [40 40 80];
wfolder='';
load(wannier_filename,'-mat');
szw=size(wannierValues);
if length(szw)==2
    % Fix for 2D maps
    szw(3)=1;
end;
sizeWannier=szw(1:3);

nBands = nOrbitals*N^2
%if (~exist('N','var'))
%    N1 = sqrt(nBands/nOrbitals);
%    if ~(N1==N)
%        disp(['Error: wrong settings for N:',num2str(N),' vs. ', num2str(N1)]);
%    end
%end;
if (~exist('xrange','var'))
    % change default behavior: calculate the whole grid
    xrange=N/2
end;
if ~ischar(xrange)
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
numlm=(2*(ceil(N/2)-1)+1)^2;
%wAcc=zeros(numlm,nOrbitals);
lhalf=ceil(N/2)-1;
mhalf=ceil(N/2)-1;
    xGridRange=xGridRange_limits(1):xGridRange_limits(2);
    yGridRange=yGridRange_limits(1):yGridRange_limits(2);
[xmesh,ymesh]=meshgrid(xGridRange_limits(1):xGridRange_limits(2),yGridRange_limits(1):yGridRange_limits(2));
else
    load(xrange)
     if ~exist('xGridRange','var')
         xGridRange=xmesh;
     end;
     if ~exist('yGridRange','var')
         yGridRange=ymesh;
     end;
end;
for en=1:le
% put all energy dependend code here!
if nargin < 2
    set_ldosfilename;
%LDOSfileName0 = [casestring,'_Vimp_', num2str(Vimp),  '_N_', num2str(N)];
%LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_ita_', num2str(ita),'_e_',num2str(E)];
else
    % input of another filename also accepted
    LDOSfileName=ldosflnm
end;
% arkward workaround
E1=E;
load(LDOSfileName,'-mat')
E=E1;
if diagonal_GF
    diagonal_string='diag';
    latticeGreens=diag(diag(latticeGreens));
else
    diagonal_string='';
end;
for zGridPoint = zGridRange
    [a,b,c]=fileparts(LDOSfileName);
    if ~isempty(a)
        a=[a,filesep];
    end;
     if ~isempty(wfolder)
        wfolder=[wfolder,filesep];
    end;
    if ~(exist([a,wfolder],'dir'))
        mkdir([a,wfolder])
    end;
    ldosmapfilename=[a,wfolder,b,c,'_z_',num2str(zGridPoint),diagonal_string];
    if exist(ldosmapfilename, 'file') == 2
        disp([ldosmapfilename,' already exists, not calculating, please delete it before recalculating it.']);
    else
localLdos=0*xmesh;
    disp(['Calculating ',num2str(zGridPoint), 'of (',num2str(zGridRange(1)),'..',num2str(zGridRange(numel(zGridRange))),')']);
    for s=1:numel(xmesh)
    %localLdos = zeros(length(xGridRange),length(yGridRange));
    %countLoopX = 0;
    %for xGridPoint = xGridRange_limits(1):xGridRange_limits(2)
    xGridPoint=xmesh(s);
    if mod(s,size(xmesh,2))==0
        disp(['Calculating ',num2str(xmesh(s)), ' of ',num2str(size(xmesh,2))]);
    end;
     %   countLoopX = countLoopX + 1;
      %  countLoopY = 0;
        %for yGridPoint =yGridRange_limits(1):yGridRange_limits(2)
        yGridPoint=ymesh(s);
         %   countLoopY = countLoopY + 1;
            r = [xGridPoint, yGridPoint, zGridPoint];
           % wAcc = [];
           %position=0;
           wAcc=zeros(numlm*nOrbitals,1);
           % change way of calculation to reduce loop
           % calculation of minimal / maximal l and m where we have Wannier
           % boxes
           centerm=yGridPoint/RDiscrete(1);
           centerl=xGridPoint/RDiscrete(1);
           minl=round(centerl-(sizeWannier(1)-shift(1)+0.5)/RDiscrete(1)+0.5);
           maxl=round(centerl+(sizeWannier(1)-shift(1)+0.5)/RDiscrete(1)-0.5);
           minm=round(centerm-(sizeWannier(2)-shift(2)+0.5)/RDiscrete(2)+0.5);
           maxm=round(centerm+(sizeWannier(2)-shift(2)+0.5)/RDiscrete(2)-0.5);
           lrange=(minl:maxl);
           mrange=(minm:maxm);
            for l = lrange
                for m = mrange
            %for l = -(ceil(N/2)-1):(ceil(N/2)-1)%lrange
             %   for m = -(ceil(N/2)-1):(ceil(N/2)-1)
                    %position=position+1;
                    position=N*(mod(l+lhalf,N))+mod(m+mhalf,N)+1;
                    %if position>numlm
                    %    disp(['outside']);
                    %end;
                    %R = [l m 0];
                    %latticeVector = RDiscrete.*R;
                    %wannierArgument = r - latticeVector;
                    %shiftedArgument = wannierArgument + shift; % translate wannier origin
                    % shift to central area if needed (for x and y direction)
                    shiftedArgument = r -  RDiscrete.*[l m 0] + shift;
                    % check whether this argument is in range or not
                    %if (((shiftedArgument) <= sizeWannier) & ((shiftedArgument) >= [1 1 1]))
                        % yes in range, now find the value
                        %wannierValue = wannierI(shiftedArgument);
                        %for orbital = 1:nOrbitals
                            w = squeeze(wannierValues(shiftedArgument(1),shiftedArgument(2),shiftedArgument(3),:));
                            wAcc(((position-1)*nOrbitals+1):position*nOrbitals) =  w;
                        %end;
                    %else
                  %disp(['outside:',num2str(shiftedArgument)]);
                        % not in range, set it to zero
                       % wAcc = [wAcc; zeros(nOrbitals,1)];
                   % end
                end
            end
         %   wAcc1=wAcc(:);
            localLdos(s) = (-1/pi)*imag(wAcc'*(latticeGreens*wAcc));
        end
    end
    % output of result
    % save the output together with geometry information necessary to plot
    save(ldosmapfilename,'localLdos','xGridRange','yGridRange','shift','sizeWannier','RDiscrete','sublattice');
    end
end
end
