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
    xGridRange = round(RDiscrete(1)/2*(xrange(1)-1)):round(RDiscrete(1)/2*xrange(2));
    %yGridRange = -60:80;
    yGridRange = round(RDiscrete(2)/2*yrange(1)):round(RDiscrete(2)/2*(yrange(2)+1));
elseif sublattice==0
    xGridRange = round(RDiscrete(1)*xrange(1)):round(RDiscrete(1)*xrange(2));
    %yGridRange = -60:80;
    yGridRange = round(RDiscrete(2)*yrange(1)):round(RDiscrete(2)*yrange(2));
end;
wannierValuesreshape=reshape(wannierValues,sizeWannier(1)*sizeWannier(2)*sizeWannier(3),nOrbitals);
wannierValuesreshape=permute(wannierValuesreshape,[2,1]);
%fileName = ['supercell_local_ldos_FeSe_U_0955_Vimp_5','_N_',num2str(N),'_M_',num2str(M),'_E_',num2str(E),'_ita_',num2str(ita)];
numlm=(2*(ceil(N/2)-1)+1)^2;
%wAcc=zeros(numlm,nOrbitals);
lhalf=ceil(N/2)-1;
mhalf=ceil(N/2)-1;
for zGridPoint = zGridRange
    disp(['Calculating ',num2str(zGridPoint), 'of (',num2str(zGridRange(1)),'..',num2str(zGridRange(numel(zGridRange))),')']);
    localLdos = zeros(length(xGridRange),length(yGridRange));
    countLoopX = 0;
    for xGridPoint = xGridRange
        disp(['Calculating ',num2str(xGridPoint), 'of (',num2str(xGridRange(1)),'..',num2str(xGridRange(numel(xGridRange))),')']);
        countLoopX = countLoopX + 1;
        countLoopY = 0;
        for yGridPoint = yGridRange
            countLoopY = countLoopY + 1;
            r = [xGridPoint, yGridPoint, zGridPoint];
            % wAcc = [];
            %position=0;
            wAcc=zeros(nOrbitals,numlm);
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
            [l1,m1]=meshgrid(lrange,mrange);
            l1=l1(:);
            m1=m1(:);
            szl=numel(l1);
            position=N*(mod(l1+lhalf,N))+mod(m1+mhalf,N)+1;
            shiftr=repmat(shift,szl,1);
            rr=repmat(r,szl,1);
            RDiscreter=repmat(RDiscrete,szl,1);
            shiftedArgument = rr -  RDiscreter.*[l1 m1 zeros(szl,1)] + shiftr;
            indexshiftedArgument=shiftedArgument(:,1)+(shiftedArgument(:,2)-1)*sizeWannier(1)+(shiftedArgument(:,3)-1)*sizeWannier(1)*sizeWannier(2);
            % for n=1:szl
            % remove loop totally
            %for l = lrange
            %   for m = mrange
            %for l = -(ceil(N/2)-1):(ceil(N/2)-1)%lrange
            %   for m = -(ceil(N/2)-1):(ceil(N/2)-1)
            %position=position+1;
            %position=N*(mod(l1(n)+lhalf,N))+mod(m1(n)+mhalf,N)+1;
            %if position>numlm
            %    disp(['outside']);
            %end;
            %R = [l m 0];
            %latticeVector = RDiscrete.*R;
            %wannierArgument = r - latticeVector;
            %shiftedArgument = wannierArgument + shift; % translate wannier origin
            % shift to central area if needed (for x and y direction)
            % check whether this argument is in range or not
            %if (((shiftedArgument) <= sizeWannier) & ((shiftedArgument) >= [1 1 1]))
            % yes in range, now find the value
            %wannierValue = wannierI(shiftedArgument);
            %for orbital = 1:nOrbitals
            %                 w = squeeze(wannierValues(shiftedArgument(n,1),shiftedArgument(n,2),shiftedArgument(n,3),:));
            %  w1=squeeze(wannierValuesreshape(:,indexshiftedArgument));
            wAcc(:,position) = wannierValuesreshape(:,indexshiftedArgument);
            %end;
            %else
            %disp(['outside:',num2str(shiftedArgument)]);
            % not in range, set it to zero
            % wAcc = [wAcc; zeros(nOrbitals,1)];
            % end
            %end
            %end
            %  wAcc=wAcc';
            wAcc1=wAcc(:);
            localLdos(countLoopX,countLoopY) = (-1/pi)*imag(wAcc1'*(latticeGreens*wAcc1));
        end
    end
    % output of result
    ldosmapfilename=[LDOSfileName,'_z_',num2str(zGridPoint),diagonal_string];
    % save the output together with geometry information necessary to plot
    save(ldosmapfilename,'localLdos','xGridRange','yGridRange','shift','sizeWannier','RDiscrete','sublattice');
end
