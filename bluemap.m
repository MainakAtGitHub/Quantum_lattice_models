function p=bluemap(figure1,offset)
if nargin < 2
    offset=0;
end;
%color1=[1 0 0]; % red
%color2=[0 0 0]; % black
color2=[1 1 1]; % white
color3=[0 0 1]; % blue
x=offset/(1-offset);
y=1/(1-offset);
input1=([-10 0 1 10]+x)/y;
colormatrix=[color2;color2;color3; color3];
maxabsekkn=1;
x=0:maxabsekkn/255:maxabsekkn;
r=interp1(input1*maxabsekkn,colormatrix(:,1),x);
g=interp1(input1*maxabsekkn,colormatrix(:,2),x);
b=interp1(input1*maxabsekkn,colormatrix(:,3),x);
p=[r',g',b'];
set(figure1,'Colormap', p);
