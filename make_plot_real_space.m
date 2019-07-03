function make_plot_real_space(inputfile)

Displaynames={'$d_{xy}$','$d_{x^2-y^2}$','$d_{xz}$','$d_{yz}$','$d_{x^2-y^2}$'};

     read_input_file=inputfile;
     read_input;
     read_input_file

if ~exist('sublattice','var')
    sublattice=input('Sublattice: ');
end;
plot_real_space(BdGfileName,'s',Displaynames,nOrbitals,sublattice)