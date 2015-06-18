function figure1=plot_wannier(wannierfile,isovalue)
load(wannierfile);
[x,y,z]=meshgrid(xGrid,yGrid,zGrid);
wsize=size(wannierValues,4);
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
xlabel({'x [bohr]'});

% Create ylabel
ylabel({'y [bohr]'});

% Create zlabel
zlabel({'z [bohr]'});

% Create textbox
annotation(figure1,'textbox',...
    [0.584984375 0.809895833333333 0.295875 0.0989583333333333],...
    'String',{['Isovalue ',num2str(isovalue),' bohr^{-3/2}']},...
    'FitBoxToText','off',...
    'LineStyle','none');
 print('-djpeg', ['/tmp/Wannier',num2str(isovalue),'band_',num2str(band),'.jpg'],'-r200');
 % plot also a cut in 2D
 figure2=figure;
 z0=34;
 surf(x(:,:,z0),y(:,:,z0),wannierValues(:,:,z0,band),'LineStyle','none');
 mx=max(max(abs(wannierValues(:,:,64,band))));
 caxis([-mx,mx]);
 blue_red_map(figure2)
axis equal
%view([0 90]);
% bug in matlab (does not show 2D objects with view([0 90]);
view([0 90.1]);
% Create xlabel
xlabel({'x [bohr]'});

% Create ylabel
ylabel({'y [bohr]'});
annotation(figure2,'textbox',...
    [0.656357142857143 0.123809523809524 0.147214285714286 0.0738095238095241],...
    'String',{['z= ',num2str(z(1,1,z0)),' bohr']},...
    'FitBoxToText','off','LineStyle','none');
% plot line around elementary cell
xline=[shift(1),shift(1)+RDiscrete(1),shift(1)+RDiscrete(1),shift(1),shift(1)]-RDiscrete(1)/2;
yline=[shift(2),shift(2),shift(2)+RDiscrete(2),shift(2)+RDiscrete(2),shift(2)]-RDiscrete(2)/2;
hold on
plot3(xGrid(xline),yGrid(yline),ones(5,1)*2*mx,'k','LineStyle','--','LineWidth',2);

 print('-djpeg', ['/tmp/Wannier_band_',num2str(band),'_z_',num2str(z0),'.jpg'],'-r200');

end;