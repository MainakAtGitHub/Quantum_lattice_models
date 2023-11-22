function f=plot_bandstructure_BZ2(inputfile,k_points)


     read_input_file=inputfile;
     read_input;
     read_input_file

% set up a mesh in the BZ
[kx,ky]=meshgrid([-pi:pi/k_points:pi]);
% put the k-points into the correct vectorlike format
klist=[kx(:),ky(:)];
% calculate the eigenvalues in the normal state for all k points
[kSpaceHamiltonian]=bands_precalc_klist2(inputfile,klist);
% put back the result into a matrix form
% To do: loop over the number of orbitals
% think about how the number of orbitals are related to the size of the
% variable kSpaceHamiltonian
f=figure;
% for orb=1:nOrb
% here some change is needed for multi orbital systems
for it_orb=1:nOrbitals
%     f=figure;
    Hk=reshape(kSpaceHamiltonian(:,it_orb),2*k_points+1,2*k_points+1);
    surf(kx,ky,Hk,it_orb*4*ones(size(Hk)));
    hold on
end
xlabel('kx');
ylabel('ky');
zlabel('band energy');
% calculating band structure along high symmetry path below
clrlst=['r','g','b','c','m']

GBZ = [0,0];
XBZ=[pi,0];
MBZ=[pi,pi];
k_points=5*k_points;
klistGXx=GBZ(1):pi/k_points:XBZ(1);
klistGXy=zeros(size(klistGXx));
klistXMy=XBZ(2):pi/k_points:MBZ(2);
klistXMx=pi*ones(size(klistXMy));
klistMGx=MBZ(1):-pi/k_points:GBZ(1);
klistMGy=MBZ(2):-pi/k_points:GBZ(2);

klist1=[klistGXx,klistXMx,klistMGx;klistGXy,klistXMy,klistMGy]';

sympath=zeros(1,length(klist1(:,1)));
sympath(1)=0;%sqrt((klist1(1,1))^2+(klist1(1,2))^2);
for i=2:length(klist1(:,1))
    sympath(i)=sympath(i-1) + sqrt((klist1(i,1)-klist1(i-1,1))^2+(klist1(i,2)-klist1(i-1,2))^2);
end

[kSpaceHamiltonian]=bands_precalc_klist2(inputfile,klist1);
f = figure;

% Create axes
axes1 = axes('Parent',f);
hold(axes1,'on');

for it_orb=1:nOrbitals
%     plot(sympath,kSpaceHamiltonian(:,it_orb),['.',clrlst(it_orb)]);hold on;%plot(sympath,kSpaceHamiltonian(:,it_orb),'*');
      plot(sympath,kSpaceHamiltonian(:,it_orb),'.','Color',[1,0,0]);hold on;
end
set(axes1,'FontSize',30,'XTick',[0 pi 2*pi 2*pi+pi*sqrt(2)],'XTickLabel',...
    {'\Gamma','X','M','\Gamma'});
xlim([0,(2+sqrt(2))*pi]);
xbnd=xlim;ybnd=ylim;
% hold on;plot(xbnd,[0,0],'k');%plot([0,(2+sqrt(2))*pi],[0,0],'k');
% hold on;plot([2*pi,2*pi],ybnd,'k');%plot([2*pi,2*pi],[-0.3,0.3],'k');
% hold on;plot([pi,pi],ybnd,'k');%plot([pi,pi],[-0.3,0.3],'k');

% xlabel('high sym path (\Gamma to X to M to \Gamma)');
ylabel('E(k) (eV)')
title(['n=',num2str(n0),', kT=',num2str(kT),' eV']);
testgca=gca;
testgca.FontSize=30;
xlim(xbnd);ylim(ybnd);
axis square;
box on; grid on;





% end

% To do: set labels of the plot with the xlabel ylabel and zlabel command!