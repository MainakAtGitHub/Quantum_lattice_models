function [ wannierValues xGrid yGrid zGrid] = process_wannier( filename, nOrb)
readstring=['%g %g %g '];
for n=1:nOrb
    readstring=[readstring,' (%g,%g)'];
end;
fid = fopen(filename);
WF = fscanf(fid, readstring, [3+2*nOrb inf]);
fclose(fid);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
szwf=size(WF);
num=0;
while ~(num==szwf(2))
    xpoints=input('points in x-direction: ');
    ypoints=input('points in y-direction: ');
    zpoints=input('points in z-direction: ');
    num=xpoints*ypoints*zpoints;
end
for n=1:nOrb
   % wannierValues(:,n)=WF(3+(n-1)*2+1,:)+1i*WF(3+n*2,:);
   % ignore the complex part
   wannierValues(:,n)=WF(3+(n-1)*2+1,:);
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
wv=real(wannierValues(:,:,:,1));
wannierValues=wv;
shift = [51 51 51];
%sizeWannier = [101 101 81];
RDiscrete = [20 20 200];
save([filename,'_conv.mat'],'wannierValues','xGrid','yGrid','zGrid','RDiscrete','shift');
end

