function f=plot_bandstructure_BZ(inputfile,k_points)

% set up a mesh in the BZ
[kx,ky]=meshgrid([-pi:pi/k_points:pi]);
% put the k-points into the correct vectorlike format
klist=[kx(:),ky(:)];
% calculate the eigenvalues in the normal state for all k points
[kSpaceHamiltonian]=bands_precalc_klist(inputfile,klist);
% put back the result into a matrix form
% To do: loop over the number of orbitals
% think about how the number of orbitals are related to the size of the
% variable kSpaceHamiltonian
f=figure
% for orb=1:nOrb
% here some change is needed for multi orbital systems
Hk=reshape(kSpaceHamiltonian,2*k_points+1,2*k_points+1);
surf(kx,ky,Hk);
hold on
% end

% To do: set labels of the plot with the xlabel ylabel and zlabel command!