function g=generate_egrid(step,maxe)
if nargin <1
step=0.05;
end;
if nargin <2
    maxe=0.03;
end;
E=[-(1.+step/2):step:(1+step/2)]*maxe
save(['egrid_step_',num2str(step),'_maxe_',num2str(maxe),'.mat'],'E');
