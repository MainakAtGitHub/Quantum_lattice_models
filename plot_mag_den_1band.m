function [ output_args ] = plot_mag_den_1band( inputfile )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
        read_input_file=inputfile;
        read_input;
        read_input_file
        load(pos_file,'-mat');
        load(BdGfileName);
        [filepath,name,ext]=fileparts(BdGfileName);
        
        figure;scatter(r(:,1),r(:,2),15+200*(nUp-nDown).^2,(nUp-nDown),'filled');
        colorbar;
        axis equal;
        title('Magnetization');
        print_pdf([filepath,filesep,name,'_magnetization.pdf']);
        
        figure;scatter(r(:,1),r(:,2),15+30*(nUp+nDown).^2,(nUp+nDown),'filled');
        colorbar;
        axis equal;
        title('Density');
        print_pdf([filepath,filesep,name,'_density.pdf']);
end

