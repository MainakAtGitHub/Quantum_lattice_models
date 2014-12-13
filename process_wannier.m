function [ wannierValues xGrid yGrid zGrid] = process_wannier( filename, nOrb,im)
readstring=['%g %g %g '];
% only two dimensional map
%readstring=['%g %g '];
if nargin <3
    im=false
end;
if im
    imagnum=2;
else
    imagnum=1;
end;
for n=1:nOrb
    if im
        readstring=[readstring,' (%g,%g)'];
        % new format
       % readstring=[readstring,' (%g,%g)'];
    else
        readstring=[readstring,' %g'];
    end;
end;


fid = fopen(filename);
% 3D maps
pos=3;
% 2D maps
%pos=2;
dble=1;
%WF = fscanf(fid, readstring, [3+imagnum*nOrb inf]);
% new format (2D)
WF = fscanf(fid, readstring, [pos+imagnum*nOrb*dble inf]);
fclose(fid);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
szwf=size(WF);
num=0;
while ~(num==szwf(2))
    disp(['total points ',num2str(szwf(2))]);
    xpoints=input('points in x-direction: ');
    ypoints=input('points in y-direction: ');
    zpoints=input('points in z-direction: ');
    num=xpoints*ypoints*zpoints;
end
shift = [51 51 51];
a=input('shift in x-direction: ');
if ~isempty(a)
    shift(1)=a;
    
end;
a=input('shift in y-direction: ');
if ~isempty(a)
    shift(2)=a;
end;
a=input('shift in z-direction: ');
if ~isempty(a)
    shift(3)=a;
end;
RDiscrete = [40 40 80];
a=input('cell size in x-direction: ');
if ~isempty(a)
    RDiscrete(1)=a;
end;
a=input('cell size in y-direction: ');
if ~isempty(a)
    RDiscrete(2)=a;
end;
a=input('cell size in z-direction: ');
if ~isempty(a)
    RDiscrete(3)=a;
end;
for n=1:nOrb
   % wannierValues(:,n)=WF(3+(n-1)*2+1,:)+1i*WF(3+n*2,:);
   % ignore the complex part
   if im
       wannierValues(:,n)=WF(pos+(n-1)*2+1,:);
   else
       wannierValues(:,n)=WF(pos+(n-1)+1,:);
   end;
end;
wannierValues=reshape(wannierValues,xpoints,ypoints,zpoints,nOrb);
minx=min(WF(1,:));
miny=min(WF(2,:));
minz=min(WF(3,:));
maxx=max(WF(1,:));
maxy=max(WF(2,:));
maxz=max(WF(3,:));
xGrid=minx:(maxx-minx)/(xpoints-1):maxx;
yGrid=miny:(maxy-miny)/(ypoints-1):maxy;
zGrid=minz:(maxz-minz)/(zpoints-1):maxz;
% some code for the BSCCO input
if nOrb==1
wv=real(wannierValues(:,:,:,1));
wannierValues=wv;
end;
%wannierValuest=wannierValues;
% for i=1:nOrb
%     for z=1:zpoints
%         wannierValues(:,:,z,i)=wannierValuest(:,:,z,mod(i+nOrb/2-1,nOrb)+1);
%     end;
% end;
%sizeWannier = [101 101 81];
%RDiscrete = [40 40 80];
save([filename,'_conv_a.mat'],'wannierValues','xGrid','yGrid','zGrid','RDiscrete','shift');
end

