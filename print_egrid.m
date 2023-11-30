function outfile=print_egrid(inputfile)
% derived from ldos_map.m r47
% inputs
%nOrbitals = 10;
%M = 5;
%ita = .001;
% default for FeSe
sublattice=1;
read_input_file=inputfile;
read_input;
read_input_file
if ~(exist('Greensenergy','var'))
    disp('No Greensenergy given, setting to 0.');
    E=0;
else
    if ~ischar(Greensenergy)
        E=Greensenergy;
    else
        % load the list of energies to be calculated from the given
        % file
        load(Greensenergy,'E');
        le=numel(E);
    end;
end
 outfile=tempname;
diary(outfile)
diary on
for n=1:le
    disp(num2str(E(n)));
end;
diary off
disp(outfile);

