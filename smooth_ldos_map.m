function s=smooth_ldos_map(inputfile,smoothparameter)
if nargin < 2
    smoothparameter=1;
end;
read_input_file=inputfile;
read_input;
read_input_file;
E = Greensenergy;
set_ldosfilename;
% some default values for the wannier functions
shift = [51 51 41];
sizeWannier = [101 101 81];
RDiscrete = [40 40 80];
if (~exist('wannier_filename','var'))
    wannier_filename='wannier_FeSe_4d_matrix_v2.mat'
end;
load(wannier_filename,'-mat');
if (~exist('diagonal_GF','var'))
    diagonal_GF=false;
end;
if exist('loacalLdos','var')
    localLdos = loacalLdos;
    clear loacalLdos;
end;
if diagonal_GF
    diagonal_string='diag';
else
    diagonal_string='';
end;
ldosmapfilename=[LDOSfileName,'_z_',num2str(0),diagonal_string];
load(ldosmapfilename,'-mat');
if exist('loacalLdos','var')
    localLdos = loacalLdos;
    clear loacalLdos;
end;
szldos=size(localLdos);
localldos3=zeros(szldos(1),szldos(2),sizeWannier(3));
localldos3(:,:,1)=localLdos;
% read in all ldos maps from z=1 to z=sizeWannier(3)
for z=2:sizeWannier(3)
        disp(['Reading in ',num2str(z), 'of ',num2str(sizeWannier(3))]);
            ldosmapfilename=[LDOSfileName,'_z_',num2str(z-(sizeWannier(3)+1)/2),diagonal_string];
    load(ldosmapfilename,'-mat');
    if exist('loacalLdos','var')
        localLdos = loacalLdos;
        clear loacalLdos;
    end;
    localldos3(:,:,z)=localLdos;
end;
% smooth data
        disp(['smoothing...']);
%[localldos3s,s] = smoothn(localldos3,smoothparameter);
localldos3s=smooth_g(localldos3,smoothparameter,[0,0,0]);
% write out data
for z=1:sizeWannier(3)
    disp(['Writing out ',num2str(z), 'of ',num2str(sizeWannier(3))]);
    ldosmapfilename=[LDOSfileName,'_z_',num2str(z-(sizeWannier(3)+1)/2),diagonal_string,'smooth'];
    localLdos=localldos3s(:,:,z);
    save(ldosmapfilename,'localLdos','xGridRange','yGridRange');
end;