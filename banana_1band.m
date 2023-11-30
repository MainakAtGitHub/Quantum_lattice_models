function [kxmin,kymin]=banana_1band(inputfile,energy)
% calculate the positions of the ends of the banana shaped equienergylines
% for a 1-band SC (assuming that the result of the homogeneous calculation
% has been written out)

% load the input file to set the variables, gave up the old .mat file
% format
read_input_file=inputfile;
read_input;
read_input_file
M=N*M;
infile=[BdGfileName,'_M_', num2str(M)];
%save([outfile1,'_normal'], 'kSpaceEigenValuesNormal', 'kx', 'ky', '-mat');
load([infile,'_SC'],'-mat');
% cut out first quadrant to speed up calculation and avoid multiples of pi
% in the result
kx=kx(1:numel(kx)/2);
ky=ky(1:numel(ky)/2);
kSpaceEigenValues=kSpaceEigenValues(1:numel(kx),1:numel(kx),:);
kxmin=inf*ones(numel(energy),1);
kymin=kxmin;
for n=1:numel(energy)
    % calculate the contour for the given energies
    if energy(n) >0
        cntr=contourc(kx,ky,kSpaceEigenValues(:,:,2),[energy(n) energy(n)]);
    else
        cntr=contourc(kx,ky,kSpaceEigenValues(:,:,1),[energy(n) energy(n)]);
    end
    % find the smallest value of k_y in the result
    % note: the output cntr is a set as follows
    % cnt=[ energy(n) numel1 el1x el1y el2x el2y ... elnumel1x elnumel1y
    %       energy(n) numel2 el1x el1y el2x el2y ... elnumel1x elnumel1y
    %       ....]
    cntr=cntr(:);
    m=3;
    while m< numel(cntr)
        % get the k-values
        kvalues=cntr(m:m+cntr(m-1)*2-1);
        % separate kx and ky
        kval=reshape(kvalues,2,numel(kvalues)/2);
        kymin1=min(kval(2,:));
        if kymin1 < kymin(n)
            kymin(n)=kymin1;
            kxmin(n)=kval(1,find(kval(2,:)==kymin1,1));
        end;
        m=m+cntr(m-1)*2+2;
    end;
end;

