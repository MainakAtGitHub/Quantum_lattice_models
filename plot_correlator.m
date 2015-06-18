function f=plot_correlator(inputfile,scale,tickx,nOrbitals,sublattice)
if nargin < 1
    inputfile='BdG_homogeneous_FeSe_Toms_BS_6Dec13_GammaCut_3_N_9.mat'
end;
if nargin <2
    scale = 's';
end;
if nargin < 3 
    %tickx={'$d_{z^2}$','$d_{x^2-y^2}$','$d_{yz}$','$d_{xz}$','$d_{xy}$'}; % orbitals order for FeSe (Tom)
% one band case
tickx={''}
end;
if nargin <4
    nOrbitals=1;
end;
if nargin <5
    sublattice=0;
end;
if nOrbitals==1
orbitalstring=''
else
orbitalstring='^{\mu\nu}'
end;
fsz=20;
% load the lattice Greens function
load(inputfile,'-mat'); % BdG_homogeneous_FeSe_Toms_BS_6Dec13_GammaCut_3_N_9.mat
N=sqrt(size(latticeGreens,1)/nOrbitals);
ldos=-1/pi*imag(latticeGreens);
for i=1:N
    for j=1:N
        diagldos(i,j,:,:)=ldos((i-1)*N*nOrbitals+(j-1)*nOrbitals+1:(i-1)*N*nOrbitals+(j-1)*nOrbitals+nOrbitals,(i-1)*N*nOrbitals+(j-1)*nOrbitals+1:(i-1)*N*nOrbitals+(j-1)*nOrbitals+nOrbitals);
    end
end;
%diagldos=diag(ldos);
for a=1:nOrbitals
    for b=1:nOrbitals
        diff((1:N)+(a-1)*N,(1:N)+(b-1)*N)=diagldos(:,:,a,b)-mean(mean(diagldos(:,:,a,b)));
    end
end;
%diagldos=reshape(diagldos,N*nOrbitals,N*nOrbitals);
%diff=(diagldos1-mean(diagldos1(:)));
cptn=['G_{RR}',orbitalstring,' [eV]^{-1}'];
% test plot to show how homogeneous the solution is
[~,h,~]=realspaceplot(diff,N,tickx,'mean1234',scale,cptn);
center=N*N/2+0.5;
for i=1:N
    for j=1:N
        offdiagdos(i,j,:,:)=latticeGreens((i-1)*N*nOrbitals+(j-1)*nOrbitals+1:(i-1)*N*nOrbitals+(j-1)*nOrbitals+nOrbitals,(1:nOrbitals)+center*nOrbitals);
    end
end;
%diagldos=diag(ldos);
for a=1:nOrbitals
    for b=1:nOrbitals
                offdiagdos1((1:N)+(a-1)*N,(1:N)+(b-1)*N)=offdiagdos(:,:,a,b);
    end
end;
%offdiagdos=latticeGreens(:,1225);
% we assume an odd number of sites in x and y directions
%offdiagGF=reshape(latticeGreens(:,center),N,N);
offdiagGF=offdiagdos1;
cptn=['Re G_{R_0 R}',orbitalstring,'[eV]^{-1}'];
[~,h,~]=realspaceplot(real(offdiagGF),N,tickx,[inputfile,'real1234'],scale,cptn);
cptn=['Im G_{R_0 R}',orbitalstring,'[eV]^{-1}'];
[~,h,~]=realspaceplot(imag(offdiagGF),N,tickx,[inputfile,'imag1234'],scale,cptn);