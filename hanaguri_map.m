function p=hanaguri_map(figure1,offset)
if nargin < 2
    offset=0;
end;
%color1=[1 0 0]; % red
color0=[ 0 0 0];
color1=[ 0 0 0]; % black
color2=[ 52 64 145]/255; % blue
color3=[137 194 66]/255; % green
color4=[239 233 60]/255; % yellow
color5=[197 38 41]/255; % red
color6=[1  1 1]; % white
x=offset/(1-offset);
y=1/(1-offset);
input1=([-10 0 0.25 0.5 0.75 1 10]+x)/y;
colormatrix=[color0; color1;color2; color3;color4;color5; color6];
maxabsekkn=1;
x=0:maxabsekkn/255:maxabsekkn;
r=interp1(input1*maxabsekkn,colormatrix(:,1),x);
g=interp1(input1*maxabsekkn,colormatrix(:,2),x);
b=interp1(input1*maxabsekkn,colormatrix(:,3),x);
set(figure1,'Colormap', [r',g',b']);
