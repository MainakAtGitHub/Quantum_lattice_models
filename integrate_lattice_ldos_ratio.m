function i=integrate_lattice_ldos_ratio(inputfile,emax,step,emin)
% integration of the ldos maps from emin to emax with step
% simple routine that reads ldos maps and adds them together
if nargin < 3
    step=1
end;
if nargin < 2
    emax=10
end
if nargin < 4
    emin=0
end
read_input_file=inputfile;
read_input;
ldossummed=[];
    read_input_file;
for i=emin:step:emax
    % read the input
    % set the energy
    E = i/1000;
    set_ldosfilename;
    ldosfilename=[LDOSfileName];
    load(ldosfilename,'-mat');
    latticeGreensp=latticeGreens;
    E =- i/1000;
    set_ldosfilename;
    ldosfilename=[LDOSfileName];
    load(ldosfilename,'-mat');    
    % add ldos maps together
    ldosfilename=[LDOSfileName,'_rel',num2str(-E)];
    latticeGreens=latticeGreensp./latticeGreens;
    E=-E;
    save(ldosfilename,'latticeGreens','N','nOrbitals','E');
    if isempty(ldossummed)
        ldossummed=latticeGreens;
    else
        ldossummed=ldossummed+latticeGreens;
    end;
end;
latticeGreens=ldossummed;
% write output
%shift(3)=-1;
ldosfilename=[LDOSfileName,'int_rel',num2str(emax)];
save(ldosfilename,'latticeGreens','N','nOrbitals');

