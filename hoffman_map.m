function p=hoffman_map(figure1,offset)
if nargin < 2
    offset=0;
end;
%color1=[1 0 0]; % red
color0=[ 1 1 1];
color1=[ 252 235 170]/255; % black
color2=[ 243 173 74]/255;
color3=[ 114 70 167]/255; % yellow
color4=[61 75 167]/255; % bright yellow
color5=[122 211 220]/255; % bright yellow
color6=[16 14 12]/255; % bright yellow

x=offset/(1-offset);
y=1/(1-offset);
input1=([-10 0 0.15 0.3 0.45 0.60 0.75 10]+x)/y;
colormatrix=[color0; color1;color2; color3;color4; color5; color6; color6];
maxabsekkn=1;
x=0:maxabsekkn/255:maxabsekkn;
r=interp1(input1*maxabsekkn,colormatrix(:,1),x);
g=interp1(input1*maxabsekkn,colormatrix(:,2),x);
b=interp1(input1*maxabsekkn,colormatrix(:,3),x);
set(figure1,'Colormap', [r',g',b']);
