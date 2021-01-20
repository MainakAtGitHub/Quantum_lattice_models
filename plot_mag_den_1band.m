function [ output_args ] = plot_mag_den_1band( inputfile,orbno )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

        if nargin < 2
            orbno=1;
        end

        read_input_file=inputfile;
        read_input;
        read_input_file
        load(pos_file,'-mat');
        load(BdGfileName);
        [filepath,name,ext]=fileparts(BdGfileName);
        
        if ~exist('nAnoUpDown','var')
            nAnoUpDown=zeros(size(nUp));
            nAnoDownUp=zeros(size(nUp));
        end
        if ~exist('field','var')
            field=[0,0,0];
        end
        
        figure;
        blue_red_map(gcf);
        set(gcf, 'Position', get(0, 'Screensize'));
        scatter(r(:,1),r(:,2),0*10+500*(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)).^2,(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)),'filled');
        hold on;scatter(r(:,1),r(:,2),[],(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)),'LineWidth',1);
   
%         if exist('nAnoUpDown','var')
            hold on;
            scl=1;
            Max_abs_mag = max(sqrt(sum([(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)).^2,(1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).^2,(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)).^2],2)))
            Avg_abs_magnetization = sum(sqrt((nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)).^2+(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)).^2+(1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).^2))/numel(nAnoUpDown(orbno:nOrbitals:end))
            Avg_vec_mag = sum([(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)),(1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))),(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end))],1)/numel(nAnoUpDown(orbno:nOrbitals:end))
%             Avg_vec_xy_mag = sum([(nAnoUpDown+nAnoDownUp).^2,(1i*(nAnoUpDown-nAnoDownUp)).^2,(nUp-nDown).^2],1)/numel(nAnoUpDown)
            Avg_xy_abs_magnetization = sum(sqrt((nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)).^2+(1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).^2))/numel(nAnoUpDown(orbno:nOrbitals:end))
            Max_xy_abs_mag = max(sqrt(sum([(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)).^2,(1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).^2],2)))           
            quiver(r(:,1),r(:,2),(scl/Max_xy_abs_mag)*(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)),(scl/Max_xy_abs_mag)*1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end)),0,'color','black');            
%         end
        
        colorbar;
%         c=colorbar;c.Label.String='M_{z}';
        axis equal;
        box on;
%       axis off;
        xlim([-N/2,N/2]);
        ylim([-N/2,N/2]);
        
%       title([name,'\n','magnetization','\n']);

        if ~exist('nAnoUpDown','var')
            title({' ',[repmat(' ',1,20),'Magnetization',' orbital/layer ',num2str(orbno)],[repmat(' ',1,20),' (',' kT=',num2str(kT),' U=',num2str(U),' n=',num2str(n0),' N=',num2str(N),' B=',num2str(field),' Mz=',num2str(sum(nUp-nDown)/numel(nUp)),' )'],' '});
        else
            title({' ',[repmat(' ',1,20),'Magnetization',' orbital/layer ',num2str(orbno),' (',' kT=',num2str(kT),' U=',num2str(U),' n=',num2str(n0),' N=',num2str(N),' B=',num2str(field),')'],[' |M|xy max=',num2str(Max_xy_abs_mag),', |M| max=',num2str(Max_abs_mag)],['Avg M=','(',num2str(Avg_vec_mag(1)),',',num2str(Avg_vec_mag(2)),',',num2str(Avg_vec_mag(3)),')',', Avg |M|=',num2str(Avg_abs_magnetization),],' '});
        end
        testgca=gca;
        testgca.FontSize=20;%30;
        axis equal;
        box on;
        xlim([-N/2,N/2]);
        ylim([-N/2,N/2]);

        print_pdf([filepath,filesep,name,'_orbital_',num2str(orbno),'_magnetization.pdf']);
        
        figure;
        scatter(r(:,1),r(:,2),[],(((nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end))).^2+(1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).^2).^0.5,'filled');
        
        figure;
        blue_red_map(gcf);
        set(gcf, 'Position', get(0, 'Screensize'));
        quiver3(r(:,1),r(:,2),zeros(size(r(:,1))),(scl/Max_abs_mag)*(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)),(scl/Max_abs_mag)*1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end)),(scl/Max_abs_mag)*(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)),0,'color','black');
        colorbar;
        testgca=gca;
        testgca.FontSize=30;%20;
        axis equal;
        box on;
        xlim([-N/2,N/2]);
        ylim([-N/2,N/2]);
        print_pdf([filepath,filesep,name,'_orbital_',num2str(orbno),'_3d_magnetization.pdf']);

        
        
%         if exist('nAnoUpDown','var')
            abs_magnetization = sqrt((nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)).^2+(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)).^2+(1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).^2);
            figure;
            bluemap(gcf);

            set(gcf, 'Position', get(0, 'Screensize'));
            scatter(r(:,1),r(:,2),1000*abs_magnetization.^2,abs_magnetization,'filled');
            hold on; scatter(r(:,1),r(:,2),[],abs_magnetization);
            colorbar;
            axis equal;
            box on;
% %             xlim([-10,10]);
% %             ylim([-10,10]);

            xlim([-N/2,N/2]);
            ylim([-N/2,N/2]);
% %             axis off;
% %             ax = gca;
% %             ax.XAxisLocation = 'origin';
% %             ax.YAxisLocation = 'origin';

            title({' ',[repmat(' ',1,20),' Orbital/layer ',num2str(orbno),' absolute magnetization (',' kT=',num2str(kT),' U=',num2str(U),' n=',num2str(n0),' N=',num2str(N),' )'],' '});
            testgca=gca;
            testgca.FontSize=20;%30;
            print_pdf([filepath,filesep,name,'_orbital_',num2str(orbno),'_abs_mag.pdf']);
%         end
%        orient landscape; saveas(gcf,[filepath,filesep,name,'_magnetization.pdf']);
%         if exist('nAnoDownUp','var')
%             figure;
%             %         blue_red_map(gcf);
%             set(gcf, 'Position', get(0, 'Screensize'));
%             quiver(r(:,1),r(:,2),(nAnoUpDown+nAnoDownUp),1i*(-nAnoUpDown+nAnoDownUp));
%             hold on;%scatter(r(:,1),r(:,2),[],(nUp-nDown),'LineWidth',1);
% 
%             colorbar;
%             axis equal;
%             %       title([name,'\n','magnetization','\n']);
%             title({' ',[repmat(' ',1,20),'MagnetizationXY (',' kT=',num2str(kT),' U=',num2str(U),' n=',num2str(n0),' N=',num2str(N),' )'],' '});
%             print_pdf([filepath,filesep,name,'_magnetizationXY.pdf']);
%             %        orient landscape; saveas(gcf,[filepath,filesep,name,'_magnetization.pdf']);
%         end

        
        figure;
        bluemap(gcf);
        set(gcf, 'Position', get(0, 'Screensize'));
        scatter(r(:,1),r(:,2),0*5+4^2*(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)).^2,(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)),'filled');
        hold on;scatter(r(:,1),r(:,2),[],(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)),'LineWidth',1);
        
        colorbar;
        axis equal;
        box on;
        xlim([-N/2,N/2]);
        ylim([-N/2,N/2]);

%       title([name,'\n','density','\n']);
        %title({' ',[repmat(' ',1,20),'Orbital/layer ',num2str(orbno),' density (',' kT=',num2str(kT),' U=',num2str(U),' n=',num2str(n0),' N=',num2str(N),' )'],' '});
        testgca=gca;
        testgca.FontSize=20;%30;
        print_pdf([filepath,filesep,name,'_orbital_',num2str(orbno),'_density.pdf']);
%        orient landscape; saveas(gcf,[filepath,filesep,name,'_density.pdf']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%combined mag den below
        figure;
%         blue_red_map(gcf);
        set(gcf, 'Position', get(0, 'Screensize'));
        
        s(1)=subplot(1,2,1);
        
        
        scatter(r(:,1),r(:,2),0*5+500*(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)).^2,(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)),'filled');
        hold on;scatter(r(:,1),r(:,2),[],(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)),'LineWidth',1);
        
        colorbar;
        axis equal;
        box on;
        xlim([-N/2,N/2]);
        ylim([-N/2,N/2]);

%       title([name,'\n','magnetization','\n']);
        %title(s(1),{' ',[repmat(' ',1,20),'Orbital/layer ',num2str(orbno),' magnetization (',' kT=',num2str(kT),' U=',num2str(U),' n=',num2str(n0),' N=',num2str(N),' )'],' '});
        
        s(2)=subplot(1,2,2);
        
        
        scatter(r(:,1),r(:,2),0*5+4^2*(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)).^2,(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)),'filled');
        hold on;scatter(r(:,1),r(:,2),[],(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)),'LineWidth',1);
        
        colorbar;
        axis equal;
        box on;
        xlim([-N/2,N/2]);
        ylim([-N/2,N/2]);

%       title([name,'\n','density','\n']);
        %title(s(2),{' ',[repmat(' ',1,20),'Orbital/layer ',num2str(orbno),' density (',' kT=',num2str(kT),' U=',num2str(U),' n=',num2str(n0),' N=',num2str(N),' )'],' '});
        
        print_pdf([filepath,filesep,name,'_orbital_',num2str(orbno),'_magnetization,density.pdf']);
%        orient landscape; saveas(gcf,[filepath,filesep,name,'_magnetization,density.pdf']);

%%%%%%%%%%%%%%%%%%%%%%%%%
       figure;
       set(gcf, 'Position', get(0, 'Screensize'));
       
       s(1)=subplot(1,3,1);
%        h1=subplot(1,3,1)
%        s(1)=blue_red_map(gca);       
       scatter(r(:,1),r(:,2),0*10+0.5^2*500*(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)).^2,(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)),'filled');
       hold on;scatter(r(:,1),r(:,2),10*ones(size(r(:,1))),(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)),'LineWidth',0.1);
       hold on;
       quiver(r(:,1),r(:,2),(scl/Max_xy_abs_mag)*(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)),(scl/Max_xy_abs_mag)*1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end)),0,'color','black');
       cb1=colorbar;
       axis equal;
       box on;
       xlim([-N/2,N/2]);
       ylim([-N/2,N/2]);
       testgca=gca;
       testgca.FontSize=10;
       axis equal;
       box on;
       xlim([-N/2,N/2]);
       ylim([-N/2,N/2]);
       title('(a)','Position', [0.25, -22, 0]);
       
       s(2)=subplot(1,3,2);
%        h2=subplot(1,3,2)
%        s(2)=bluemap(gca);
       scatter(r(:,1),r(:,2),0.5^2*1000*abs_magnetization.^2,abs_magnetization,'filled');
       hold on;scatter(r(:,1),r(:,2),10*ones(size(r(:,1))),abs_magnetization,'LineWidth',0.1);
       cb2=colorbar;
       axis equal;
       box on;
       xlim([-N/2,N/2]);
       ylim([-N/2,N/2]);
       testgca=gca;
       testgca.FontSize=10;
       title('(b)','Position', [0.25, -22, 0]);
       
       s(3)=subplot(1,3,3);
%        h3=subplot(1,3,3)
%        s(3)=bluemap(gca);
       scatter(r(:,1),r(:,2),0*5+0.5^2*4^2*(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)).^2,(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)),'filled');
       hold on;scatter(r(:,1),r(:,2),10*ones(size(r(:,1))),(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)),'LineWidth',0.1);       
       cb3=colorbar;
       axis equal;
       box on;
       xlim([-N/2,N/2]);
       ylim([-N/2,N/2]);
       testgca=gca;
       testgca.FontSize=10;
       title('(c)','Position', [0.25, -22, 0]);
       
       print_pdf([filepath,filesep,name,'_orbital_',num2str(orbno),'_mag_abs_mag_density.pdf']);




end


