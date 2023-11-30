function p=blue_red_map(figure1,offset)
if nargin < 2
    offset=0;
end
color1=[1 0 0]; % red
%color2=[0 0 0]; % black
color2=[1 1 1]; % white
color3=[0 0 1]; % blue
x=real(offset)/(1-real(offset));
y=1/(1-real(offset));
input1=([-10 -1 0 1 10]+x)/y;
colormatrix=[color1; color1; color2;color3; color3];
maxabsekkn=1;
x=-maxabsekkn:maxabsekkn/127:maxabsekkn;
if imag(offset)>0
    x=x(end:-1:1);
end
r=interp1(input1*maxabsekkn,colormatrix(:,1),x);
g=interp1(input1*maxabsekkn,colormatrix(:,2),x);
b=interp1(input1*maxabsekkn,colormatrix(:,3),x);
p=[r',g',b'];
if nargin>0
    set(figure1,'Colormap', [r',g',b']);
end
