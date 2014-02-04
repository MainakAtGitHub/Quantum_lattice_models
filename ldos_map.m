function r=ldos_map(inputfile,ldosflnm)

% inputs
%nOrbitals = 10;
%M = 5;
%ita = .001;
read_input_file=inputfile;
read_input;
read_input_file
if (~exist('xrange','var'))
    xrange=3
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
N = sqrt(nBands/nOrbitals);

% Local greens function
% Wannier vector
yrange=xrange;
%xGridRange = -80:60;
xGridRange = -RDiscrete(1)/2*(xrange+1):RDiscrete(1)/2*xrange;
%yGridRange = -60:80;
yGridRange = -RDiscrete(1)/2*yrange:RDiscrete(1)/2*(yrange+1);

%fileName = ['supercell_local_ldos_FeSe_U_0955_Vimp_5','_N_',num2str(N),'_M_',num2str(M),'_E_',num2str(E),'_ita_',num2str(ita)];
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
            wAcc = []; 
            for l = -(ceil(N/2)-1):(ceil(N/2)-1)
                for m = -(ceil(N/2)-1):(ceil(N/2)-1)
                    R = [l m 0];
                    latticeVector = RDiscrete.*R;
                    wannierArgument = r - latticeVector;
                    shiftedArgument = wannierArgument + shift; % translate wannier origin
                    % check whether this argument is in range or not
                    if (((shiftedArgument) <= sizeWannier) & ((shiftedArgument) >= [1 1 1]))
                        % yes in range, now find the value
                        %wannierValue = wannierI(shiftedArgument);
                        %for orbital = 1:nOrbitals
                            w = squeeze(wannierValues(shiftedArgument(1),shiftedArgument(2),shiftedArgument(3),:));
                            wAcc = [wAcc; w];
                        %end;
                    else
                        % not in range, set it to zero
                        wAcc = [wAcc; zeros(nOrbitals,1)];
                    end
                end
            end
            localLdos(countLoopX,countLoopY) = (-1/pi)*imag(wAcc'*(latticeGreens*wAcc));
        end
    end
    % output of result
    ldosmapfilename=[LDOSfileName,'_z_',num2str(zGridPoint),diagonal_string];
    % save the output together with geometry information necessary to plot
    save(ldosmapfilename,'localLdos','xGridRange','yGridRange','shift','sizeWannier','RDiscrete');
end
