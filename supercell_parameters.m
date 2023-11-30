function [BZ]= supercell_parameters(M)
% set the k-grid and the weights for a supercell calculation
[kx,ky]=meshgrid(0:(M-1));
kx=kx(:)*(2*pi/M);
ky=ky(:)*(2*pi/M);
BZ.k=[kx,ky];
BZ.weight=ones(1,M^2)/M^2;
end
