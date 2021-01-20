function [kSpaceHamiltonian]=Hk_klist1(inputfile,klist)

% for given inputfile (*.txt) and a file containing k-points, calculate the momentum space Hamiltonian
% of the band structure

% by default set the chemical potential to zero
 mu=0;
    % load the input file to set the variables, gave up the old .mat file
    % format
     read_input_file=inputfile;
     % call the script (not function that parses the input file and sets up
     % the correct variables)
     read_input;
    % read_input_file


% load the hoppings from the file
% load(TB_file,'-mat');
% nOrbitals = size(TBparameters,1);
latticeVector=[1,0,0;
                -1 0 0;
                0 1 0;
                0 -1 0;
                1 1 0;
                -1 -1 0;
                1 -1 0;
                -1 1 0;
                2 0 0;
                -2 0 0;
                0 2 0;
                0 -2 0;
                0 0 0];
if bilayer_int
    perpTB=0.08*0.15*[0;0;0;0;-1/2;-1/2;-1/2;-1/2;1/4;1/4;1/4;1/4;1];
end
TBparameters=zeros(nOrbitals,nOrbitals,size(latticeVector,1));
if exist('ref_grid_hopping_file','var')
    t_r_ref=load(ref_grid_hopping_file);
end
for it_latticeVector=1:size(latticeVector,1) % hopping 0 for the last vector [0 0 0], hence the -1
%     for it_orb=1:nOrbitals
        if it_latticeVector~=size(latticeVector,1)
            TBparameters(1,1,it_latticeVector) = griddata(t_r_ref(:,1),t_r_ref(:,2),(10^(-3))*t_r_ref(:,3),...
                latticeVector(it_latticeVector,1),latticeVector(it_latticeVector,2),'natural');
            TBparameters(2,2,it_latticeVector)=TBparameters(1,1,it_latticeVector);
        end
        if bilayer_int
            TBparameters(1,2,it_latticeVector)=perpTB(it_latticeVector);
            TBparameters(2,1,it_latticeVector)=perpTB(it_latticeVector);
%             TBparameters(1,2,it_latticeVector)=0.08*TBparameters(1,1,it_latticeVector);
%             TBparameters(2,1,it_latticeVector)=0.08*TBparameters(1,1,it_latticeVector);
        end
%     end
end
% if bilayer_int
%     TBparameters(1,2,end) = 0.08*0.15; %t_perp; NN hopping = 0.15 eV;
%     TBparameters(2,1,end) = 0.08*0.15;
% end

% number of hoppings in band structure
nUnitCells = size(latticeVector,1);

% allow to set the chemical potential
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) = ...
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) - mu*eye(nOrbitals);

% fix for the lattice vectors (no 3rd component allowed)
latticeVector=latticeVector(:,1:2);

% read in k values if klist is a string, else just use the klist as data
if isa(klist,'char') || isa(klist,'string') 
	kin=load(klist);
else
    kin=klist;
end
szk=size(kin);
% initialize the result variables
kSpaceHamiltonian = zeros(szk(1), 2*nOrbitals, 2*nOrbitals);

%iterate over the momentum vectors
for iKx = 1:szk(1)
	    % set the momentum vector (in 2D)
            k = [kin(iKx,1) kin(iKx,2)];
	    % initialize the Hamiltonian to zeros
            kSpaceHopping = zeros(nOrbitals,nOrbitals);
	    % add all terms coming from each hopping
            for iUnitCell = 1:nUnitCells
                iLatticeVector = latticeVector(iUnitCell,:);
		% do (manually) the Fourier transform and add to the Hamiltonian
                kSpaceHopping = kSpaceHopping + TBparameters(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
                if bilayer_int
%                                   kSpaceHopping(1,2) = 0.08*0.15; %t_perp; NN hopping = 0.15 eV;
%                                   kSpaceHopping(2,1) = 0.08*0.15;
%                                   kSpaceHopping(1,2) = 0.08*0.15*(cos(k(1))-cos(k(2)))^2;
%                                   kSpaceHopping(2,1) = 0.08*0.15*(cos(k(1))-cos(k(2)))^2;
                end
            end
            % put in the result into the multidimensional array
            if int_soc
                term_soc_1=[0.03*0.15*(sin(k(2)+1i*sin(k(1)))),0;0,-0.03*0.15*(sin(k(2)+1i*sin(k(1))))];
                term_soc_2=[0.03*0.15*(sin(k(2)-1i*sin(k(1)))),0;0,-0.03*0.15*(sin(k(2)-1i*sin(k(1))))];
            else 
                term_soc_1=zeros(size(kSpaceHopping));
                term_soc_2=zeros(size(kSpaceHopping));
            end
            kSpaceHamiltonian(iKx, :,:) =  [kSpaceHopping,term_soc_1;
                                           term_soc_2,kSpaceHopping]; %lambda_soc=0.03*0.15; 
end
