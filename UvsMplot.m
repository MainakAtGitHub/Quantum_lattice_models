function f=UvsMplot(folderName)
% BdGfolder='/home/UFAD/mainak.pal/Desktop/summer_2019_office/andreas/BdG';
BdGfolder = '/home/UFAD/mainak.pal/Desktop/summer_2019_office/andreas/BdG/Updated_synced_BdG_scripts_21June2020';
cd(folderName);
inFileListStruc=dir('*.txt');
plotList=[];
cd(BdGfolder);
for i=1:length(inFileListStruc)
    itrFile=inFileListStruc(i).name;
    inputfile=[num2str(folderName),'/',num2str(itrFile)];
    read_input_file=inputfile;
    read_input;
    read_input_file
    load(BdGfileName);
    absMagNorm=sum(abs(nUp-nDown))/length(nUp);
    signedMagNorm=sum((nUp-nDown))/length(nUp);
    plotList=[plotList;[U,absMagNorm,signedMagNorm]];
end
[srtdx,srtdxI]=sort(plotList(:,1));
plotList(:,1)=plotList(srtdxI,1);
plotList(:,2)=plotList(srtdxI,2);
save([folderName,'/','N_',num2str(N),'kT_',num2str(kT),'n_',num2str(n0),'_plotList.mat'],'plotList');
% save([folderName,'/','N_',num2str(N),'kT_',num2str(kT),'n_',num2str(n0),'_plotList.txt'],'plotList','-ascii');
figure;

% scatter(plotList(:,1),plotList(:,2),50+plotList(:,2).^2,plotList(:,2),'b','filled');

% scatter(plotList(:,1),plotList(:,2),100+plotList(:,2).^2,[0.6 0.3 0],'filled');
plot(plotList(:,1),plotList(:,2),'LineWidth',2,'color',[0.5 0.5 0]);hold on;
plot(plotList(:,1),plotList(:,2),'.','MarkerSize',30,'color',[0.5 0.5 0]);
% legend('Homogeneous/RandomIniGuess/size30/n0\_0.95','Location','Best');
title({'','U vs avg. absolute M',''});
xlabel('U');
ylabel('avg. absolute M');
xlim([0.33-0.0001,0.6+0.0001]);ylim([-.0001,0.3]);grid off;
testgca=gca;
testgca.FontSize=20;%20;
daspect([1,1.6,2]);
% axis equal;
box on;
% axis equal;
% colorbar;
set(gcf, 'Position', get(0, 'Screensize'));
% grid on;
saveas(gcf,[folderName,'/','N_',num2str(N),'kT_',num2str(kT),'n_',num2str(n0),'_UvsAvgAbsM.fig']);
saveas(gcf,[folderName,'/','N_',num2str(N),'kT_',num2str(kT),'n_',num2str(n0),'_UvsAvgAbsM.jpg']);
orient landscape;
saveas(gcf,[folderName,'/','N_',num2str(N),'kT_',num2str(kT),'n_',num2str(n0),'_UvsAvgAbsM.pdf']);


% hold on;

% % % % % figure;
% % % % % 
% % % % % % scatter(plotList(:,1),plotList(:,2),50+plotList(:,2).^2,plotList(:,2),'b','filled');
% % % % % scatter(plotList(:,1),plotList(:,3),100+plotList(:,2).^2,[0.3 0.6 0],'filled');
% % % % % % legend('Homogeneous/RandomIniGuess/size30/n0\_0.95','Location','Best');
% % % % % title({'','U vs avg. signed M',''});
% % % % % xlabel('U');
% % % % % ylabel('avg. signed M');
% % % % % % axis equal;
% % % % % % colorbar;
% % % % % set(gcf, 'Position', get(0, 'Screensize'));
% % % % % grid on;
% % % % % saveas(gcf,[folderName,'/','N_',num2str(N),'kT_',num2str(kT),'n_',num2str(n0),'_UvsSignedAbsM.fig']);
% % % % % saveas(gcf,[folderName,'/','N_',num2str(N),'kT_',num2str(kT),'n_',num2str(n0),'_UvsSignedAbsM.jpg']);
% % % % % orient landscape;
% % % % % saveas(gcf,[folderName,'/','N_',num2str(N),'kT_',num2str(kT),'n_',num2str(n0),'_UvsSignedAbsM.pdf'])
% % % % % 
% % % % % % hold on;
end