function r=filter_ldos_map(inputfile,ldosflnm)

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
if (~exist('Greensenergy','var'))
    Greensenergy=0.0084
end;
if (~exist('STMpositionfile','var'))
    % defines file that contains the postions of the TIP to be calculated
    % as well as the average area + weights
    STMpositionfile='positions.mat'
end;
if (~exist('diagonal_GF','var'))
    diagonal_GF=false;
end;
if diagonal_GF
    diagonal_string='diag';
    latticeGreens=diag(diag(latticeGreens));
else
    diagonal_string='';
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

shift = [51 51 41];
RDiscrete = [40 40 80];
load(STMpositionfile,'-mat');
%szw=size(wannierValues);
%sizeWannier=szw(1:3);

%nBands = size(latticeGreens,1);

if (~exist('xrange','var'))
    % change default behavior: calculate the whole grid
    xrange=N/2
end;

% be a little more general to allow for maps that are not centered at the
% central site
if numel(xrange)==1
    xrange(2)=xrange;
    xrange(1)=-xrange(1);
end;
if (~exist('yrange','var'))
    yrange=xrange;
end;
% size of the weight
szw=size(stmweight);
% number of positons
npos=size(position,1);
% initialize the result
STM_ldos=zeros(1,npos);
for zGridPoint = zGridRange
        ldosmapfilename=[LDOSfileName,'_z_',num2str(zGridPoint),diagonal_string];
    load(ldosmapfilename,'-mat');
    disp(['Integrating map locally for z= ',my_int2str(zGridPoint)]);
    szlldos=size(localLdos);
    for p=1:size(position,1)
        % the coordinate system for the position is fixed to the center of
        % the current map
        % calculate the center of the STM tip in terms of discrete indices
        discreteposition=round(position(p,:).*RDiscrete(1:2));
        % cut out the localLdos that is contributing
        posmin=round(szlldos/2)-fix(szw/2)+discreteposition;
        posmax=round(szlldos/2)+fix(szw/2)+discreteposition;
        localLdos_cut=localLdos(posmin(1):posmax(1),posmin(2):posmax(2));
    % take the convolution of the map with the average area
    STM_ldos(p) =sum(sum( localLdos_cut.*stmweight));
    end
    save([ldosmapfilename,'weight'],'STM_ldos');
end
