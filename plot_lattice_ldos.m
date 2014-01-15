% input
load lattice_greens_supercell_FeSe_N_15_M_9_U_0955_Vimp_5_E_minPt0084.mat
N = 15;
E = .0084;
nOrbitals = 10;


% total LDOS at Fe sites
latticeGreensDiag = diag(latticeGreens);
ldos = (-1/pi)*imag(reshape(latticeGreensDiag,nOrbitals*N,N));
Fe1LDOS = zeros(N,N);
Fe2LDOS = zeros(N,N);
for i = 1:N
    Fe1LDOS(i,:) = sum(ldos((i-1)*nOrbitals + (1:nOrbitals/2),:),1);
    Fe2LDOS(i,:) = sum(ldos((i-1)*nOrbitals + ((nOrbitals/2 + 1):nOrbitals),:),1);
end
Fe1LDOS = Fe1LDOS';
Fe2LDOS = Fe2LDOS';


% rotating coordinates by pi/4
n = ceil(N/2);
ldos2plot = diag(Fe1LDOS(:,n),0);
for i=1:(n-1)
    ldos2plot = ldos2plot + diag(Fe1LDOS((i+1):(N-i),n+i), 2*i) + diag(Fe1LDOS((i+1):(N-i),n-i), -2*i) + ....
    diag(Fe2LDOS(i:(N-i),n+i), 2*i-1) + diag(Fe2LDOS(i:(N-i),n-i+1), -(2*i-1));
end


% plotting
ldos2plotAppended = flipud(ldos2plot);
ldos2plotAppended(:,N+1) = 0;
ldos2plotAppended(N+1,:) = 0;
figure;
pcolor(ldos2plotAppended);
axis('square');
title(['E = ',num2str(E), ' eV']);
xlabel('x')
ylabel('y')
colorbar