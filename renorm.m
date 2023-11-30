function r=renorm(E)
E_0=2
Z_0=1/3
r=E.*(1-(1-Z_0)*exp(-E.^2/(2*E_0^2)));