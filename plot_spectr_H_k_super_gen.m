function outpt = plot_spectr_H_k_super_gen(inputfile,klist,klist_no_points)
if nargin < 3
    klist_no_points = '';
end
    [pthstr,nm,xt] = fileparts(inputfile);
    load([pthstr,'/sup_cell_k_spectr_',num2str(klist_no_points),'y','.mat'],'klist','pf_or_eig_val_list');
    fg=figure;%('Position', [400,550,500,300]);
    plot(klist(:,2)/pi, pf_or_eig_val_list,'Color','r','Marker','.' );
    xlabel('k_y/\pi');
    ylabel('Energy (eV)')
    testgca=gca;
    testgca.FontSize=15;
    grid on;
    print_pdf([pthstr,'/sup_cell_spectr_',num2str(klist_no_points),'_gap_ky.pdf']);
    saveas(gcf,[pthstr,'/sup_cell_spectr_',num2str(klist_no_points),'_gap_ky.fig']);
    
    fg1=figure;%('Position', [400,550,500,300]);
    plot(klist(:,2)/pi, pf_or_eig_val_list,'Color','r','Marker','.' );
    xlabel('k_y/\pi');
    ylabel('Energy (eV)')
    testgca=gca;
    testgca.FontSize=15;
    ylim([-0.05,0.05]);
    grid on;
    print_pdf([pthstr,'/sup_cell_spectr_',num2str(klist_no_points),'_gap_ky_zoomed.pdf']);
    
    fg2=figure;%('Position', [400,550,500,300]);
    plot(klist(:,2)/pi, pf_or_eig_val_list,'Color','r','Marker','.' );
    xlabel('k_y/\pi');
    ylabel('Energy (eV)')
    testgca=gca;
    testgca.FontSize=15;
    ylim([-0.1,0.1]);
    grid on;
    print_pdf([pthstr,'/sup_cell_spectr_',num2str(klist_no_points),'_gap_ky_zoomed_intermediate.pdf']);