function p=song_map(figure1,offset)
if nargin < 2
    offset=0;
end;
%color1=[1 0 0]; % red
color0=[62 16 16]/255;
color1=[62 16 16]/255; % black
color2=[164 52 52]/255;
color3=[241 184 94]/255;; % yellow
color4=[246 225 178]/255;; % bright yellow
color5=[1 1 1]; % white
x=offset/(1-offset);
y=1/(1-offset);
input1=([-10 0 0.33 0.66 1 10]+x)/y;
colormatrix=[color0; color1;color2; color3;color5;color5];
maxabsekkn=1;
x=0:maxabsekkn/255:maxabsekkn;
r=interp1(input1*maxabsekkn,colormatrix(:,1),x);
g=interp1(input1*maxabsekkn,colormatrix(:,2),x);
b=interp1(input1*maxabsekkn,colormatrix(:,3),x);
set(figure1,'Colormap', [r',g',b']);
