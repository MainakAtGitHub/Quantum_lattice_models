prefix='~/itp/docs/susc/numerics/cuda/nocuda/'
gapfile=[prefix,'FeSe_BdG/FeSe_tom_10_true_2Dklist_kz_0_full_BZ_large.dat1gap_FT_summed.dat'];
kfile=[prefix,'FeSe_BdG/klist_kz_0_full_BZ_large.dat'];
gapfile=[prefix,'FeSe_10band/FeSe_tom_10_true_2Dklist_kz_0_full_BZ.dat1gap_FT_summed.dat']
%kfile=[prefix,'klist_kz_0_full_BZ.dat'];
gapfile=[prefix,'FeSe_10band/FeSe_tom_10_true_2Dklist_kz_0_full_BZ_large.dat1gap_FT_summed.dat']


gap=load(gapfile);
k=load(kfile);
szk=size(k);
scale='s';
% convert to meV
gap_s=1000*reshape(gap,sqrt(szk(1)),sqrt(szk(1)));
kx=reshape(k(:,1),sqrt(szk(1)),sqrt(szk(1)))/pi;
ky=reshape(k(:,2),sqrt(szk(1)),sqrt(szk(1)))/pi;
datarealmax=max(abs(gap_s(:)));
lm=log(datarealmax)/log(10);
mtix=10^(ceil(lm));
% do some refinement to avoid only single labels
if (ceil(lm)-lm > 0.5)
    tx=[-.5:.05:.5];
elseif    (ceil(lm)-lm > 0.2)
    tx=[-1:0.1:1];
else
    tx=[-1:0.2:1];
end;
ticks=mtix*tx; 
labels = num2str(repmat(sign(ticks).*(abs(ticks)), 1, 1)', 2);
figure1=figure('Position',[100, 100, 400, 420]);
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1]);
view(axes1,[-0.5 90]);
grid(axes1,'on');
hold(axes1,'all');
if scale=='s'
    surf(kx,ky,sign(gap_s).*sqrt(abs(gap_s)),'LineStyle','none',...
    'FaceColor','interp',...
    'EdgeColor','none');
    ticks=sign(tx).*sqrt(mtix*abs(tx));
    datarealmax=sqrt(datarealmax);
else
    surf(kx,ky,gap_s,'LineStyle','none',...
    'FaceColor','interp',...
    'EdgeColor','none');
end;
ca=caxis;
caxis(([-1 1])*max(abs(ca)));
xlim([-1 1]);
ylim([-1 1]);
%xlabel('k_x/\pi');
%ylabel('k_y/\pi');
colorbar_rwb(figure1,datarealmax,ticks,labels,true);
if scale=='s'
    print('-dpng',['/tmp/test_k_gap_sqrt.png']);
else
    print('-dpng',['/tmp/test_k_gap.png']);
end;
