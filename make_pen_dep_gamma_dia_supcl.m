function [tg,superLatticeVectors] = make_pen_dep_gamma_dia_supcl(N, TBparameters, latticeVector, pen_dep_dir) % pen_dep_dir can take +,-1, +,-2
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
pen_dep_gamma_param = zeros(size(TBparameters));
for i_LV = 1: size(latticeVector,1)    
     pen_dep_gamma_param(:,:,i_LV)= -TBparameters(:,:,i_LV)*(latticeVector(i_LV,abs(pen_dep_dir))*sign(pen_dep_dir))^2;  
end
[tg,superLatticeVectors] = supercell_hoppings(N, pen_dep_gamma_param, latticeVector);
