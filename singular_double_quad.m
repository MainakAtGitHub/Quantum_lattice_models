function I = singular_double_quad(B)
% Integrates nearly singular data with integrand of type int[(1/B(kx,ky));p1,q1;p2,q2]
% INPUT:       B: 2D matrix; denominator of integrand - can contain one or more than one near
%              singular points
%           
% OUTPUT:   Integral(1/B) 
% v1.1 3/26/2012
% Ref: singular_quad

N = size(B);
B_1D = zeros(N(2),1)*1i;
    parfor j = 1:N(2)
        B_1D(j)= singular_quad(B(:,j));% integrate over first dimension
    end
I = singular_quad(1./B_1D);% integrate over second dimension
