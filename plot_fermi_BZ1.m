function [f,contours]=plot_fermi_BZ1(inputfile,k_points)

     read_input_file=inputfile;
     read_input;
     read_input_file

% set up a mesh in the BZ
kgrid=[-pi:pi/k_points:pi];
[kx,ky]=meshgrid(kgrid);
% put the k-points into the correct vectorlike format
klist=[kx(:),ky(:)];
% calculate the eigenvalues in the normal state for all k points
[kSpaceHamiltonian]=bands_precalc_klist1(inputfile,klist);
% put back the result into a matrix form
% To do: loop over the number of orbitals
% think about how the number of orbitals are related to the size of the
% variable kSpaceHamiltonian
% % % % f= figure('Position',[150, 100, 500, 500]);;

% Create axes
axes1 = axes('Parent',f);
hold(axes1,'on');

% for orb=1:nOrb
% here some change is needed for multi orbital systems
the_colors={'red','blue','green','cyan'};
for it_orb=1:2*nOrbitals
%     f=figure;
    Hk=reshape(kSpaceHamiltonian(:,it_orb),2*k_points+1,2*k_points+1);
    % plot contour lines for the Fermi surface
    contour(kgrid/pi,kgrid/pi,Hk,[0,0],'Color',the_colors{it_orb},'LineWidth',2);
    % look up the documentation of this command which essentially calculates
    % the Fermi points and understand the format of the output variable!
% %     contours=contourc(kgrid/pi,kgrid/pi,Hk,[0,0]);
    % think about this command, look up documentation
    axis equal
    hold on
end
xlabel('k_x/\pi');
ylabel('k_y/\pi');
title('fermi surface');



% set(gcf, 'Position', get(0, 'Screensize'));
% testgca=gca;
% testgca.FontSize=20;%20;

box(axes1,'on');
axis(axes1,'tight');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','DataAspectRatio',[1 1 1],'Layer','top','YTick',...
    [-1 -0.5 0 0.5 1],'YTickLabel',{'-1','-0.5','0','0.5','1'});
% end

% To do: set labels of the plot with the xlabel ylabel and zlabel command!