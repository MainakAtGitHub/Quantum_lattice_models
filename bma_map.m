function p=bma_map(figure1,offset)
if nargin < 2
    offset=0;
end;
%color1=[1 0 0]; % red
color0=[ 0 0 0];
color1=[ 0 0 1]; % black
color2=[ 0 1 0];
color3=[ 1 0 0 ]; % yellow
color4=[1  1 1]; % bright yellow
x=offset/(1-offset);
y=1/(1-offset);
input1=([-10 0 0.5 1 10]+x)/y;
colormatrix=[color0; color1;color2; color3;color4];
maxabsekkn=1;
x=0:maxabsekkn/255:maxabsekkn;
r=interp1(input1*maxabsekkn,colormatrix(:,1),x);
g=interp1(input1*maxabsekkn,colormatrix(:,2),x);
b=interp1(input1*maxabsekkn,colormatrix(:,3),x);
set(figure1,'Colormap', [r',g',b']);
