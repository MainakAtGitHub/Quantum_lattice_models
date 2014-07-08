function p=hoffman_map(figure1,offset)
if nargin < 2
    offset=0;
end;
%color1=[1 0 0]; % red
color0=[ 1 22 85]/255;
color1=[ 18 90 24]/255; % black
color2=[ 125 244 0]/255;
color3=[ 253 241 7]/255; % yellow
color4=[255 13 9]/255; % bright yellow
color5=[255 255 255]/255; % bright yellow
%color6=[16 14 12]/255; % bright yellow

x=offset/(1-offset);
y=1/(1-offset);
input1=([-10 0 0.2 0.4 0.6 0.80 1 10]+x)/y;
colormatrix=[color0; color0;color1; color2;color3; color4; color5; color5];
maxabsekkn=1;
x=0:maxabsekkn/255:maxabsekkn;
r=interp1(input1*maxabsekkn,colormatrix(:,1),x);
g=interp1(input1*maxabsekkn,colormatrix(:,2),x);
b=interp1(input1*maxabsekkn,colormatrix(:,3),x);
set(figure1,'Colormap', [r',g',b']);
