function [figure1,psi2]=plot_wannier(wannierfile,isovalue,z0)
load(wannierfile,'-mat');
wannierValues=permute(wannierValues,[2,1,3,4]);
conv=true;
if conv
zGrid=zGrid*0.529177249
xGrid=xGrid*0.529177249;
yGrid=yGrid*0.529177249;
end
[x,y,z]=meshgrid(xGrid,yGrid,zGrid);
wsize=size(wannierValues,4);
if wsize==10
    wsize=10%
end;
for band=1:wsize

fv = isosurface(x,y,z,wannierValues(:,:,:,band),isovalue);
fvm = isosurface(x,y,z,wannierValues(:,:,:,band),-isovalue);
figure1=figure;
p=patch(fv);
isonormals(x,y,z,wannierValues(:,:,:,band),p)
set(p,'FaceColor','blue','EdgeColor','none');
hold on
pm=patch(fvm);
isonormals(x,y,z,wannierValues(:,:,:,band),pm)
set(pm,'FaceColor','red','EdgeColor','none');
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud

% Create xlabel
xlabel({'x (Å)'});

% Create ylabel
ylabel({'y (Å)'});

% Create zlabel
zlabel({'z (Å)'});

if false
% Create textbox
annotation(figure1,'textbox',...
    [0.584984375 0.809895833333333 0.295875 0.0989583333333333],...
    'String',{['Isovalue ',num2str(isovalue),' bohr^{-3/2}']},...
    'FitBoxToText','off',...
    'LineStyle','none');
end
 print('-djpeg', ['/tmp/Wannier',num2str(isovalue),'band_',num2str(band),'.jpg'],'-r200');
 % plot also a cut in 2D
 figure2=figure;
 if nargin < 3
 z0=(size(wannierValues,3)-1)-30;
 z0=95;
  z0=100;
    z0=116;
 end;

  %z0=28
%imagesc(xGrid,yGrid,wannierValues(:,:,z0,band));
%set(gca,'ydir','normal')
 %surf(x(:,:,z0),y(:,:,z0),wannierValues(:,:,z0,band),'LineStyle','none');
 imagesc(x(:,1,z0),y(1,:,z0),wannierValues(:,:,z0,band));
 % calculate "normalization on plane"
 psi2=sum(sum(wannierValues(:,:,z0,band).*conj(wannierValues(:,:,z0,band))));
 mx=max(max(abs(wannierValues(:,:,z0,band))))
 caxis([-mx,mx]);
colorbar;
blue_red_map(figure2);
axis equal
view([0 90]);
% bug in matlab (does not show 2D objects with view([0 90]);
%view([0 90.1]);
% Create xlabel
xlabel({'x (Å)'});

% Create ylabel
ylabel({'y (Å)'});
annotation(figure2,'textbox',...
    [0.656357142857143 0.123809523809524 0.147214285714286 0.0738095238095241],...
    'String',{['z= ',num2str(z(1,1,z0)),' Å']},...
    'FitBoxToText','off','LineStyle','none');
% plot line around elementary cell
xline=[shift(1),shift(1)+RDiscrete(1),shift(1)+RDiscrete(1),shift(1),shift(1)]-RDiscrete(1)/2;
yline=[shift(2),shift(2),shift(2)+RDiscrete(2),shift(2)+RDiscrete(2),shift(2)]-RDiscrete(2)/2;
hold on
%plot3(xGrid(round(xline)),yGrid(round(yline)),ones(5,1)*2*mx,'k','LineStyle','--','LineWidth',2);
if wsize>1
% mark Li positions
NNNx=[0,0,0,1,1,1]*RDiscrete(1);
NNNy=[0,1,-1,0,1,-1]*RDiscrete(2);
scatter3(xGrid(round(NNNx+shift(1)-RDiscrete(1)/4)),yGrid(round(NNNy+shift(2)-RDiscrete(1)/4)),ones(numel(NNNx),1)*mx*1,'o','MarkerEdgeColor','k')
NNNy=[-1,-1,-1,0,0,0]*RDiscrete(2);
NNNx=[0,1,-1,0,1,-1]*RDiscrete(1);
scatter3(xGrid(round(NNNx+shift(1)+RDiscrete(1)/4)),yGrid(round(NNNy+shift(2)+RDiscrete(1)/4)),ones(numel(NNNx),1)*mx*1,'square','MarkerEdgeColor','k')
end;

 print('-djpeg', ['/tmp/Wannier_band_',num2str(band),'_z_',num2str(z0),'.jpg'],'-r200');

end;