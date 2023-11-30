function tb=convert_tb_PC_AK(tbfile)
% convert a matlab type tight-binding file back to the one having one line
% per hopping (Wannier90 type) ASCI writing
% load the data
load(tbfile,'-mat')
cut=0;
% initialize array
tb=[];
sztb=size(TBparameters);
% iterate over all orbitals and all lattice vectors
for mu=1:sztb(1)
    for nu=1:sztb(2)
        for r=1:sztb(3)
            if abs(TBparameters(mu,nu,r))>cut
                % add each nonzero value into a row of the tb array
                tb=[tb;latticeVector(r,1:2),0,mu,nu,TBparameters(mu,nu,r)];
            end
        end
    end
end
% write out the numeric result
dlmwrite([tbfile,'.csv'],tb,'precision',10);
