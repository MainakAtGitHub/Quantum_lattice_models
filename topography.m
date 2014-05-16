function s=topography(inputfile,ldosvalue,emax)
if nargin < 2
    ldosvalue=6e-4;
end;
if nargin < 3
    emax=10
end;
% read the inputvariables
read_input_file=inputfile;
read_input;
E = Greensenergy;
% put together string for ldos filename
set_ldosfilename;
% some default values for the wannier functions
shift = [51 51 41];
sizeWannier = [101 101 81];
RDiscrete = [40 40 80];
if (~exist('wannier_filename','var'))
    wannier_filename='wannier_FeSe_4d_matrix_v2.mat'
end;
% load wannier file to get the correct shift sizeWannier RDiscrete
load(wannier_filename,'-mat');
% some default values
if (~exist('diagonal_GF','var'))
    diagonal_GF=false;
end;
% bug from earlier version
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
% set up a 3D array for the ldos
plotz=(1:56)+41;
localldos3=zeros(szldos(1),szldos(2),numel(plotz));
localldos3(:,:,1)=localLdos;
% read in all ldos maps from z=1 to z=sizeWannier(3)
loopvar=0;
for z=plotz
    loopvar=loopvar+1;
    %sizeWannier(3)
            ldosmapfilename=[LDOSfileName,'_z_',num2str(z-(sizeWannier(3)+1)/2),diagonal_string,'int',num2str(emax),'smooth'];
            disp(['Reading in ',num2str(z), ' of ',num2str(plotz(end)), ' : ',ldosmapfilename]);
           %             ldosmapfilename=[LDOSfileName,'_z_',num2str(z-(sizeWannier(3)+1)/2),diagonal_string];

    load(ldosmapfilename,'-mat');
    if exist('loacalLdos','var')
        localLdos = loacalLdos;
        clear loacalLdos;
    end;
    localldos3(:,:,loopvar)=localLdos;
end;
szlL=size(localLdos);
% extrapolate values
extra=0;
if extra>0
for n=1:szlL(1)
    for m=1:szlL(2)
        localldos3(n,m,loopvar+1:loopvar+extra)=shiftdim(interp1(loopvar+[-1 0],squeeze(localldos3(n,m,loopvar+[-1 0])),loopvar+1:loopvar+extra,'linear','extrap'),-1);
    end
end;
end;
% calculate topography from the given data
s=figure;
axes1 = axes('Visible','off','Parent',s);
[x,y,z]=meshgrid(xGridRange,yGridRange,[plotz(1):plotz(end)+extra]);
z=z/80*10.44-8;
p=patch(isosurface(x,y,z,localldos3,ldosvalue));
isonormals(x,y,z,localldos3,p);
isocolors(x,y,z,z,p);
set(p,'FaceColor','interp','EdgeColor','none')
daspect([1 1 1]);axis tight
camlight; lighting none;
%isosurface(x,y,z,localldos3,ldosvalue)
xlabel('x');
ylabel('y');
zlabel('z');
view(-45,90);
% plot red lines
cut=7
sublattice=1
if sublattice==1
    offset(1)=-0.25*RDiscrete(1);
    offset(2)=0.25*RDiscrete(2);
else
    offset=[0 0];
end;
hold on
z1=[1 1]*max(z(:))*1.1;
        x1=cut/4*[-RDiscrete(1),RDiscrete(2)]+offset(1);
        y1=cut/4*[-RDiscrete(1),RDiscrete(2)]+offset(2);
        p1=plot3(x1,y1,z1,'r');
        
% Create multiple lines using matrix input to plot3
set(p1,'LineWidth',1);

        x1=cut/4*[RDiscrete(1),-RDiscrete(2)]+offset(1);
        y1=cut/4*[-RDiscrete(1),RDiscrete(2)]+offset(2);
        p2=plot3(x1,y1,z1,'r');
        set(p2,'LineWidth',1);
        % bug with openGL renderer!
        set(gcf, 'renderer', 'zbuffer');
        % another use:
        % print -zbuffer

