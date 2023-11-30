function I = singular_quad(B)
% Integrates nearly singular data with integrand of type int[(1/B(y)),p,q]
% INPUT:    B: denominator of integrand - can contain one or more than one near
%              singular points
%           
% OUTPUT:   Integral(1/B) lim: p to q
% v1.1 3/26/2012
% Ref: Tight binding models for amorphous systems ; M L Roth
% http://prb.aps.org/abstract/PRB/v7/i10/p4321_1


% M = round((length(B)-1)/2);
% I = 0;
%     for n = 1:M
%         a = .5*(B(2*n-1)-2*B(2*n)+B(2*n+1));
%         b = .5*(B(2*n+1) - B(2*n-1));
%         c = B(2*n);
%         if a ~= 0
%             x_plus = (-b + sqrt(b^2 - 4*a*c))/(2*a);
%             x_minus = (-b - sqrt(b^2 - 4*a*c))/(2*a);
%             den = (a.*(x_plus - x_minus));
%             num = log(((1 - x_plus).*(1 + x_minus))./((1 + x_plus).*(1 - x_minus)));
%             I_n = num/den;
%         else
%             if b ~= 0
%                 I_n = (1/b)*log((c + b)/(c - b));
%             else
%                 I_n = 2/c;
%             end
%         end
%         I = I + I_n;
%     end
M = length(B)-1;    
I = 0*1i;
 for n = 1:M
     a = B(n + 1) - B(n);
     b = B(n);
     if abs(a) > 1e-5
         I_n = (1/a)*log(1 + (a/b));
     else
         I_n = 1/b - a/(2*b^2) + a^2/(3*b^3);
     end
     I = I + I_n;
 end
 
 % vectorized code
 