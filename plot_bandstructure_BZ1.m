function f=plot_bandstructure_BZ1(inputfile,k_points)


     read_input_file=inputfile;
     read_input;
     read_input_file

% set up a mesh in the BZ
[kx,ky]=meshgrid([-pi:pi/k_points:pi]);
% put the k-points into the correct vectorlike format
klist=[kx(:),ky(:)];
% calculate the eigenvalues in the normal state for all k points
[kSpaceHamiltonian]=bands_precalc_klist1(inputfile,klist);
% put back the result into a matrix form
% To do: loop over the number of orbitals
% think about how the number of orbitals are related to the size of the
% variable kSpaceHamiltonian
f=figure;
% for orb=1:nOrb
% here some change is needed for multi orbital systems
for it_orb=1:2*nOrbitals
%     f=figure;
    Hk=reshape(kSpaceHamiltonian(:,it_orb),2*k_points+1,2*k_points+1);
    surf(kx,ky,Hk,it_orb*4*ones(size(Hk)));
    hold on
end
xlabel('kx');
ylabel('ky');
zlabel('band energy');
% end

% To do: set labels of the plot with the xlabel ylabel and zlabel command!