function [ output_args ] = plot_mag_den_1band( inputfile,orbno,folder_for_U_vs_Mabs_phase_data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
        if nargin<3
            folder_for_U_vs_Mabs_phase_data='';
        end
        
        if nargin < 2
            orbno=1;
        end

        read_input_file=inputfile;
        read_input;
        read_input_file
        if numel(N)==1
            N=[N,N];
        end
        if exist('pos_file','var')
            load(pos_file,'-mat');
        else
            [tempxx,tempyy]=meshgrid(-N(1)/2+0.5*mod(N(1)+1,2):N(1)/2-0.5*mod(N(1)+1,2),-N(2)/2+0.5*mod(N(2)+1,2):N(2)/2-0.5*mod(N(2)+1,2));
            r=[tempxx(:),tempyy(:)];
            clear tempxx tempyy;
        end
        load(BdGfileName);
        [filepath,name,ext]=fileparts(BdGfileName);
        
        if (nargin < 2)
            nUpdmy=zeros(length(nUp)/nOrbitals,1);
            nDowndmy=zeros(length(nDown)/nOrbitals,1);
            if exist('nAnoUpDown','var')
            nAnoUpDowndmy=zeros(length(nAnoUpDown)/nOrbitals,1);
            nAnoDownUpdmy=zeros(length(nAnoDownUp)/nOrbitals,1);
            end
            for idmy=1:nOrbitals
                nUpdmy=nUpdmy+nUp(idmy:nOrbitals:end);
                nDowndmy=nDowndmy+nDown(idmy:nOrbitals:end);
                if exist('nAnoUpDown','var')
                nAnoUpDowndmy=nAnoUpDowndmy+nAnoUpDown(idmy:nOrbitals:end);
                nAnoDownUpdmy=nAnoDownUpdmy+nAnoDownUp(idmy:nOrbitals:end);
                end
            end
            nOrbitals=1;
            nUp=nUpdmy;
            nDown=nDowndmy;
            if exist('nAnoUpDown','var')
            nAnoUpDown=nAnoUpDowndmy;
            nAnoDownUp=nAnoDownUpdmy;
            end
        end
        
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
            Max_abs_mag = max(sqrt(sum([(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)).*conj((nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end))),...
                (1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).*conj((1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end)))),...
                (nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)).*conj((nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)))],2)))
            Avg_abs_magnetization = sum(sqrt((nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)).*conj((nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)))+...
            (nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)).*conj((nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)))+...
                (1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).*conj((1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end)))))...
                )/numel(nAnoUpDown(orbno:nOrbitals:end))
            Avg_vec_mag = sum([(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)),(1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))),(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end))],1)/numel(nAnoUpDown(orbno:nOrbitals:end))
%             Avg_vec_xy_mag = sum([(nAnoUpDown+nAnoDownUp).^2,(1i*(nAnoUpDown-nAnoDownUp)).^2,(nUp-nDown).^2],1)/numel(nAnoUpDown)
            Avg_xy_abs_magnetization = sum(sqrt((nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)).*conj((nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)))+...
                (1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).*conj((1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end)))))...
                )/numel(nAnoUpDown(orbno:nOrbitals:end))
            Max_xy_abs_mag = max(sqrt(sum([(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)).*conj((nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end))),...
                (1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).*conj((1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))))]...
                ,2)))  
            Max_x_abs_mag = max(sqrt(sum([(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)).*conj((nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end))),...
                (0.*1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).*conj((1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))))]...
                ,2)))
            Max_y_abs_mag = max(sqrt(sum([0.*(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)).*conj((nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end))),...
                (1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).*conj((1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))))]...
                ,2)))
            Max_z_abs_mag = max(sqrt(sum([(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)).*conj((nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end))),...
                (0.*1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).*conj((1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))))]...
                ,2)))
            quiver(r(:,1),r(:,2),(scl/Max_xy_abs_mag)*(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)),(scl/Max_xy_abs_mag)*1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end)),0,'color','black');            
%         end
        
        colorbar;
%         caxis([-.1,.1]);
%         c=colorbar;c.Label.String='M_{z}';
        axis equal;
        box on;
%       axis off;
        xlim([-N(1)/2,N(1)/2]);
        ylim([-N(2)/2,N(2)/2]);
        
%       title([name,'\n','magnetization','\n']);

        if ~exist('nAnoUpDown','var')
             title({' ',[repmat(' ',1,20),'Magnetization',' orbital/layer ',num2str(orbno)],[repmat(' ',1,20),' (',' kT=',num2str(kT),' U=',num2str(U),' n=',num2str(n0),' N=','[',num2str(N(1)),',',num2str(N(2)),']',' B=',num2str(field),' Mz=',num2str(sum(nUp-nDown)/numel(nUp),10),' )'],' '});
        else
             title({' ',[repmat(' ',1,20),'Magnetization',' orbital/layer ',num2str(orbno),' (',' kT=',num2str(kT),' U=',num2str(U),' n=',num2str(n0),' N=','[',num2str(N(1)),',',num2str(N(2)),']',' B=',num2str(field),')'],[' |M|xy max=',num2str(Max_xy_abs_mag,10),', |M| max=',num2str(Max_abs_mag,10)],['|M|x max=',num2str(Max_x_abs_mag,10),', |M|y max=',num2str(Max_y_abs_mag,10),', |M|z max=',num2str(Max_z_abs_mag,10)],['Avg M=','(',num2str(Avg_vec_mag(1)),',',num2str(Avg_vec_mag(2)),',',num2str(Avg_vec_mag(3)),')',', Avg |M|=',num2str(Avg_abs_magnetization,10),],' '});
        end
        testgca=gca;
        testgca.FontSize=20;%20;%10;%30;
        axis equal;
        box on;
        xlim([-N(1)/2,N(1)/2]);
        ylim([-N(2)/2,N(2)/2]);
        ca=caxis; mx_ca=max(abs(ca)); caxis([-mx_ca,mx_ca]);
        print_pdf([filepath,filesep,name,'_orbital_',num2str(orbno),'_magnetization.pdf']);
        
        figure;
        scatter(r(:,1),r(:,2),[],(((nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end))).^2+(1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).^2).^0.5,'filled');
        
        figure;
        blue_red_map(gcf);
        set(gcf, 'Position', get(0, 'Screensize'));
        x1=(scl/Max_abs_mag)*(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end));
        x2=(scl/Max_abs_mag)*1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end));
        x3=(scl/Max_abs_mag)*(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end));
%         f_3d_mag=quiver3(r(:,1),r(:,2),zeros(size(r(:,1))),x1,x2,x3,'color','black');

        f_3d_mag=quiver3D([r(:,1),r(:,2),zeros(size(r(:,1)))],real([x1,x2,x3]));%,0,'color','black','LineWidth',1.5);%%%,([ones(size(r(:,1))),(r(:,2)+5)/8,zeros(size(r(:,1)))+0.1])
%         hold on; line([0.5,0.5],[-2.5,2.5],[0,0],'color','black');
%         hold on; line([-0.5,-0.5],[-2.5,2.5],[0,0],'color','black');
%         hold on; line([1.5,1.5],[-2.5,2.5],[0,0],'color','black');
        camlight head;
        lighting phong;
%         xlim([-2,2]);
        axis tight;
%         f_3d_mag.Color = 'red';
%         f_3d_mag.Marker = 'O';
        colorbar;
        testgca=gca;
        testgca.FontSize=20;%20;%20;
        axis equal;
        box on;
        xlim([-N(1)/2,N(1)/2]);
        ylim([-N(2)/2,N(2)/2]);
        print_pdf([filepath,filesep,name,'_orbital_',num2str(orbno),'_3d_magnetization.pdf']);

        
        
%         if exist('nAnoUpDown','var')
            abs_magnetization = sqrt((nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)).^2+(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)).^2+(1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end))).^2);
            figure;
            bluemap(gcf);
  %         ca=caxis;mx_ca=max(abs(ca));
%         
            if ~isempty(folder_for_U_vs_Mabs_phase_data)
                try
                    load([folder_for_U_vs_Mabs_phase_data,'/','U_vs_Mabs_phase_data.mat'],'-mat');
                catch
                    U_vs_Mabs_phase_data=[0,0];
                    save([folder_for_U_vs_Mabs_phase_data,'/','U_vs_Mabs_phase_data.mat'],'U_vs_Mabs_phase_data');
                end
                nrmlzd_abs_magnetization=sum(abs_magnetization)/numel(abs_magnetization);
                U_vs_Mabs_phase_data=[U_vs_Mabs_phase_data;[U,nrmlzd_abs_magnetization]];
                save([folder_for_U_vs_Mabs_phase_data,'/','U_vs_Mabs_phase_data.mat'],'U_vs_Mabs_phase_data' );
            end

            set(gcf, 'Position', get(0, 'Screensize'));
            scatter(r(:,1),r(:,2),100*abs_magnetization.^2,abs_magnetization,'filled');
            hold on; scatter(r(:,1),r(:,2),[],abs_magnetization);
            colorbar;
%             caxis([0,0.25]);
            axis equal;
            box on;
% %             xlim([-10,10]);
% %             ylim([-10,10]);

        xlim([-N(1)/2,N(1)/2]);
        ylim([-N(2)/2,N(2)/2]);
% %             axis off;
% %             ax = gca;
% %             ax.XAxisLocation = 'origin';
% %             ax.YAxisLocation = 'origin';

             title({' ',[repmat(' ',1,20),' Orbital/layer ',num2str(orbno),' absolute magnetization (',' kT=',num2str(kT),' U=',num2str(U),' n=',num2str(n0),' N=','[',num2str(N(1)),',',num2str(N(2)),']',' )'],' '});
            testgca=gca;
            testgca.FontSize=12;%20;%20;%30;
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
%         caxis([0,2.5]);
        axis equal;
        box on;
        xlim([-N(1)/2,N(1)/2]);
        ylim([-N(2)/2,N(2)/2]);

%        title([name,'\n','density','\n']);
        %title({' ',[repmat(' ',1,20),'Orbital/layer ',num2str(orbno),' density (',' kT=',num2str(kT),' U=',num2str(U),' n=',num2str(n0),' N=',num2str(N),' )'],' '});
        testgca=gca;
        testgca.FontSize=12;%20;%20;%30;
        print_pdf([filepath,filesep,name,'_orbital_',num2str(orbno),'_density.pdf']);
%        orient landscape; saveas(gcf,[filepath,filesep,name,'_density.pdf']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%combined mag den below
        figure;
%         blue_red_map(gcf);
        set(gcf, 'Position', get(0, 'Screensize'));
        
        s(1)=subplot(2,1,1);
        
        
        scatter(r(:,1),r(:,2),0*5+500*(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)).^2,(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)),'filled');
        hold on;scatter(r(:,1),r(:,2),[],(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)),'LineWidth',1);
        
        colorbar;
        axis equal;
        box on;
        xlim([-N(1)/2,N(1)/2]);
        ylim([-N(2)/2,N(2)/2]);

        ca=caxis; mx_ca=max(abs(ca)); caxis([-mx_ca,mx_ca]);
%       title([name,'\n','magnetization','\n']);
        %title(s(1),{' ',[repmat(' ',1,20),'Orbital/layer ',num2str(orbno),' magnetization (',' kT=',num2str(kT),' U=',num2str(U),' n=',num2str(n0),' N=',num2str(N),' )'],' '});
        
        s(2)=subplot(2,1,2);
        
        
        scatter(r(:,1),r(:,2),0*5+4^2*(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)).^2,(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)),'filled');
        hold on;scatter(r(:,1),r(:,2),[],(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)),'LineWidth',1);
        
        colorbar;
        axis equal;
        box on;
        xlim([-N(1)/2,N(1)/2]);
        ylim([-N(2)/2,N(2)/2]);

        ca=caxis; mx_ca=max(abs(ca)); caxis([0,mx_ca]);
%       title([name,'\n','density','\n']);
        %title(s(2),{' ',[repmat(' ',1,20),'Orbital/layer ',num2str(orbno),' density (',' kT=',num2str(kT),' U=',num2str(U),' n=',num2str(n0),' N=',num2str(N),' )'],' '});
        
        print_pdf([filepath,filesep,name,'_orbital_',num2str(orbno),'_magnetization,density.pdf']);
%        orient landscape; saveas(gcf,[filepath,filesep,name,'_magnetization,density.pdf']);

%%%%%%%%%%%%%%%%%%%%%%%%%
       figure;
       set(gcf, 'Position', get(0, 'Screensize'));
       
       s(1)=subplot(3,1,1);
%        h1=subplot(1,3,1)
%        s(1)=blue_red_map(gca);       
       scatter(r(:,1),r(:,2),0*10+0.5^2*500*(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)).^2,(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)),'filled');
       hold on;scatter(r(:,1),r(:,2),10*ones(size(r(:,1))),(nUp(orbno:nOrbitals:end)-nDown(orbno:nOrbitals:end)),'LineWidth',0.1);
       hold on;
       quiver(r(:,1),r(:,2),(scl/Max_xy_abs_mag)*(nAnoUpDown(orbno:nOrbitals:end)+nAnoDownUp(orbno:nOrbitals:end)),(scl/Max_xy_abs_mag)*1i*(nAnoUpDown(orbno:nOrbitals:end)-nAnoDownUp(orbno:nOrbitals:end)),0,'color','black');
       cb1=colorbar;
       axis equal;
       box on;
        xlim([-N(1)/2,N(1)/2]);
        ylim([-N(2)/2,N(2)/2]);
       testgca=gca;
       testgca.FontSize=12;
       axis equal;
       box on;
        xlim([-N(1)/2,N(1)/2]);
        ylim([-N(2)/2,N(2)/2]);
       title('(a)','Position', [0.25, -22, 0]);
       ca=caxis; mx_ca=max(abs(ca)); caxis([-mx_ca,mx_ca]);
       
       s(2)=subplot(3,1,2);
%        h2=subplot(1,3,2)
%        s(2)=bluemap(gca);
       scatter(r(:,1),r(:,2),0.5^2*100*abs_magnetization.^2,abs_magnetization,'filled');
       hold on;scatter(r(:,1),r(:,2),10*ones(size(r(:,1))),abs_magnetization,'LineWidth',0.1);
       cb2=colorbar;
       axis equal;
       box on;
        xlim([-N(1)/2,N(1)/2]);
        ylim([-N(2)/2,N(2)/2]);
       testgca=gca;
       testgca.FontSize=12;
       title('(b)','Position', [0.25, -22, 0]);
       ca=caxis; mx_ca=max(abs(ca)); caxis([0,mx_ca]);
       
       s(3)=subplot(3,1,3);
%        h3=subplot(1,3,3)
%        s(3)=bluemap(gca);
       scatter(r(:,1),r(:,2),0*5+0.5^2*4^2*(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)).^2,(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)),'filled');
       hold on;scatter(r(:,1),r(:,2),10*ones(size(r(:,1))),(nUp(orbno:nOrbitals:end)+nDown(orbno:nOrbitals:end)),'LineWidth',0.1);       
       cb3=colorbar;
       %caxis([5.5,6]);
       axis equal;
       box on;
        xlim([-N(1)/2,N(1)/2]);
        ylim([-N(2)/2,N(2)/2]);
       testgca=gca;
       testgca.FontSize=12;
       title('(c)','Position', [0.25, -22, 0]);
       ca=caxis; mx_ca=max(abs(ca)); caxis([0,mx_ca]);
       
       print_pdf([filepath,filesep,name,'_orbital_',num2str(orbno),'_mag_abs_mag_density.pdf']);




end


