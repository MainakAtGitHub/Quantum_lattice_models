function tb=convert_tb_PC_AK(tbfile)
load(tbfile,'-mat')
cut=0;
tb=[];
sztb=size(TBparameters);
for mu=1:sztb(1)
    for nu=1:sztb(2)
        for r=1:sztb(3)
            if abs(TBparameters(mu,nu,r))>cut
                tb=[tb;latticeVector(r,1:2),0,mu,nu,TBparameters(mu,nu,r)];
            end;
        end
    end
end
dlmwrite([tbfile,'.csv'],tb,'precision',10);
