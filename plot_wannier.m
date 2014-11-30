function figure1=plot_wannier(wannierfile,isovalue);
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
end;