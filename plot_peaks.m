if false
for vimp=0.4:.1:7
n=n+1;peaks{n}=make_plot_impurity_dos(['~/docs/real_space/BdG/calc/V_imp/LDOS_FeSe_Tom__Vimp_',num2str(vimp),'_N_13_M_10_tetra_corr'],0.005); close all;
end;
for vimp=8:1:10
n=n+1;peaks{n}=make_plot_impurity_dos(['~/docs/real_space/BdG/calc/V_imp/LDOS_FeSe_Tom__Vimp_',num2str(vimp),'_N_13_M_10_tetra_corr'],0.005); close all;
end;
vimp=[0.4:0.1:7;8:10]
figure
save('peaks_LDOS.mat','peaks','vimp');
end;

load peaks_LDOS
coloruf1=[250 	70 	22 ]/255;
coloruf2= [0 	48 	135]/255;
figure
for m=1:4
    %figure
    hold on
    switch m
        case 1 % impurity site
            clr=[0 0 0];
        case 2 % NN site
            clr=coloruf1;
        case 3
            clr=coloruf2;
        case 4
            clr=[1 1 0];
    end;
    for n=1:numel(peaks)
        peaka=peaks{n}{m};
        for s=1:numel(peaka)/2
            sz=round(peaka(2,s)*5);
            if sz >0
                plot(vimp(n),peaka(1,s),'k^','markerfacecolor',clr,'MarkerSize',sz);
            end;
        end
    end
    xlabel('V_{imp} [eV]');
    ylabel('\omega [eV]');
end;
        
        