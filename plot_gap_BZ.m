function f=plot_gap_BZ(inputfile,k_points)

% set up a mesh in the BZ
[kx,ky]=meshgrid([-pi:pi/k_points:pi]);
% put the k-points into the correct vectorlike format
klist=[kx(:),ky(:)];
% calculate the superconducting gap in momentum and orbital space
[kspacegap]=Deltak_klist(inputfile,klist);
% put back the result into a matrix form
% To do: loop over the number of orbitals (two loops here: Why?)
% think about how the number of orbitals are related to the size of the
% variable kSpaceHamiltonian
f=figure;
% for orb1=1:nOrb
% for orb2=1:nOrb
% here some change is needed for multi orbital systems
% Put the result for each orbital component into a block of the matrix and
% then do the plotting!
Delta_k=reshape(kspacegap,2*k_points+1,2*k_points+1);
surf(kx,ky,real(Delta_k));
hold on
% end
% end

% To do: also plot the imaginary part of the order parameter
% Why is this needed?

% To do: set labels of the plot with the xlabel ylabel and zlabel command!