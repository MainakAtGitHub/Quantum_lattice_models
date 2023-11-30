function [o,energy]=octett_1band(kx,ky,energy,symm)
% calculates the scattering vectors from the octett model given the kx, ky points for the ends of the
% "banananas" as calculated using "banana_1band.m"
% definition of the onset point as in Hanaguri et al. 
qx=-ky(:)/pi;
qy=-kx(:)/pi;
if nargin < 4
    symm=false;
end;
if symm
qx1=[qx;-qy;-qx;qy;-qx;qy;qx;-qy];
qy1=[qy;-qx;-qy;qx;qy;-qx;-qy;qx];
qx=qx1;
qy=qy1;
energy=[energy(:);energy(:);energy(:);energy(:);energy(:);energy(:);energy(:);energy(:)];
end;
% give back a struct with all the q-vectors

o.q1=[-qx-qx,qy-qy];
o.q2=[-qy-qx,qx-qy];
o.q3=[-qy-qx,-qx-qy];
o.q4=[-qx-qx,-qy-qy];
o.q5=[qx-qx,-qy-qy];
o.q6=[qy-qx,-qx-qy];
o.q7=[qy-qx,qx-qy];
% rotate the vectors according to the symmetry
