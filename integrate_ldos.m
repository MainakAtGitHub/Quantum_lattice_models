function i=integrate_ldos(inputfile,emax,step,emin)
% integration of the ldos maps from emin to emax with step
% simple routine that reads ldos maps and adds them together
if nargin < 3
    step=1
end;
if nargin < 2
    emax=10
end
if nargin < 4
    emin=0
end
read_input_file=inputfile;
read_input;
if (~exist('diagonal_GF','var'))
    diagonal_GF=false;
end;
if diagonal_GF
    diagonal_string='diag';
else
    diagonal_string='';
end;
ldossummed=[];
    read_input_file;
shift=[51 51 51];
    z=zGridRange(1)+shift(3);
for i=emin:step:emax
    % read the input
    % set the energy
    E = i/1000;
    set_ldosfilename;
    ldosmapfilename=[LDOSfileName,'_z_',num2str(z-shift(3)),diagonal_string]
    load(ldosmapfilename,'-mat');
    if exist('loacalLdos','var')
        localLdos = loacalLdos;
        clear loacalLdos;
    end;
    % add ldos maps together
    if isempty(ldossummed)
        ldossummed=localLdos;
    else
        ldossummed=ldossummed+localLdos;
    end;
end;
localLdos=ldossummed;
% write output
%shift(3)=-1;
ldosmapfilename=[LDOSfileName,'_z_',num2str(z-shift(3)),diagonal_string,'int',num2str(emax)];
save(ldosmapfilename,'localLdos','xGridRange','yGridRange','shift','sizeWannier','RDiscrete','sublattice');

