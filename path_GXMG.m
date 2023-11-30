function p=path_GXMG(points,filename)
% generate a path for the tetragonal system (usual path in 2D from Gamma-X-M-Gamma)
% points give the number of steps in each part of the path

% select the edge points
g=[0,0,0];
x=[0,1,0];
s=[1,1,0];
% generate equally spaced sequence of vectors
[kx1,ky1,kz1]=esp3(g,x,points);
[kx2,ky2,kz2]=esp3(x,s,points);
[kx3,ky3,kz3]=esp3(s,g,points);
% put vectors together
kx=[kx1,kx2,kx3];
ky=[ky1,ky2,ky3];
kz=[kz1,kz2,kz3];
% multiply with pi for convenient input to calculate bands
p=[kx;ky;kz]*pi;
% if filename argument is given, write out the data
if nargin>1
    csvwrite(filename,p')
end;
end

function [x,y,z]=esp3(A,B,n)
x=equalspace(A(1),B(1),n);
y=equalspace(A(2),B(2),n);
z=equalspace(A(3),B(3),n);
end

function es=equalspace(a,b,n)
if a-b==0
    es=ones(1,n)*a;
else
    es=(a:-(a-b)/n:b+(a-b)/n);
end;
end
