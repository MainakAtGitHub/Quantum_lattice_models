function [t,H_soc1,H_soc2] = hoppings(r,N_lin,figg,t_r_ref,eff_pot,further_N_cut_off,NN_cut_off,chg_map_filepath,nOrbitals,bilayer_int,int_soc)

% if nargin < 8
%     filepath
% end

if nargin < 11
    int_soc=false;
end

if nargin < 10
    bilayer_int=false;
end    


if nargin < 9
    nOrbitals=1;
end

if nargin < 2
    N_lin=0;
end

if nargin < 3
    figg=false;
end

if nargin < 4
    t_r_ref=[];
end

if nargin < 5
    eff_pot=1.5;
end

if nargin < 6
    further_N_cut_off=1.45;
end

if nargin < 7
    NN_cut_off=1.3;
end

t = zeros(size(r,1)*nOrbitals);

H_soc1 = zeros(size(r,1)*nOrbitals);
H_soc2 = zeros(size(r,1)*nOrbitals);

Q_eff_bondwise=zeros(size(r,1));
d=zeros(size(r,1));
d_x=zeros(size(r,1));
d_y=zeros(size(r,1));
Q_eff=zeros(size(r,1),1);
Q_Cu=+4;
Q_O=-2;
lmbda=1;
A=0.5*exp(1/4);

%%%%%%%%% setting up effective charge and effective potential in presence of dislocation

[nearest_neighbours,nearest_neighbours_index] = choose_near_neighbours(r,N_lin,NN_cut_off);
[further_neighbours,further_neighbours_index] = choose_near_neighbours(r,N_lin,further_N_cut_off);%,figg);

if bilayer_int
intlayerNN_cut_off=2.05;
[intlayerNN,intlayerNN_index] = choose_near_neighbours(r,N_lin,intlayerNN_cut_off);
 perpTBlatticeVector=[1,0,0;
                -1,0,0;
                0,1,0;
                0,-1,0;
                1,1,0;
                -1,-1,0;
                1,-1,0;
                -1,1,0;
                2,0,0;
                -2,0,0;
                0,2,0;
                0,-2,0;
                0,0,0];
            perpTB=0.08*0.15*[0;0;0;0;-1/2;-1/2;-1/2;-1/2;1/4;1/4;1/4;1/4;1];
            perpTBlatticeVector=[perpTBlatticeVector;
                3,0,0;
                -3,0,0;
                0,3,0;
                0,-3,0;
                2,2,0;
                2,-2,0;
                -2,2,0;
                -2,-2,0];
            perpTB=[perpTB;0;0;0;0;0;0;0;0];
end
if figg
    figure;
end
for i = 1:size(r,1)
    [min_abs_d_x,min_abs_d_x_index]=min([abs(r(i,1)-r(nearest_neighbours_index{i},1)),...
        abs(r(i,1)-r(nearest_neighbours_index{i},1)-N_lin),...
        abs(r(i,1)-r(nearest_neighbours_index{i},1)+N_lin)],[],2);
    
    dummy_x_dist_colmns=-[(r(i,1)-r(nearest_neighbours_index{i},1)),...
        (r(i,1)-r(nearest_neighbours_index{i},1)-N_lin),...
        (r(i,1)-r(nearest_neighbours_index{i},1)+N_lin)];
    
    for it_min=1:length(min_abs_d_x_index)
        min_signed_d_x(it_min)=dummy_x_dist_colmns(it_min,min_abs_d_x_index(it_min));
    end
    
    for it_neigh_index=1:numel(nearest_neighbours_index{i})
        d_x(i,nearest_neighbours_index{i}(it_neigh_index)) = min_signed_d_x(it_neigh_index);
    end
    
    [min_abs_d_y,min_abs_d_y_index]=min([abs(r(i,2)-r(nearest_neighbours_index{i},2)),...
        abs(r(i,2)-r(nearest_neighbours_index{i},2)-N_lin),...
        abs(r(i,2)-r(nearest_neighbours_index{i},2)+N_lin)],[],2);
    
    dummy_y_dist_colmns=-[(r(i,2)-r(nearest_neighbours_index{i},2)),...
        (r(i,2)-r(nearest_neighbours_index{i},2)-N_lin),...
        (r(i,2)-r(nearest_neighbours_index{i},2)+N_lin)];
    
    for it_min=1:length(min_abs_d_y_index)
        min_signed_d_y(it_min)=dummy_y_dist_colmns(it_min,min_abs_d_y_index(it_min));
    end
    
    for it_neigh_index=1:numel(nearest_neighbours_index{i})
        d_y(i,nearest_neighbours_index{i}(it_neigh_index)) = min_signed_d_y(it_neigh_index);
    end
    
        
    d(i,nearest_neighbours_index{i}) = (d_x(i,nearest_neighbours_index{i}).^2 + d_y(i,nearest_neighbours_index{i}).^2).^(0.5);
    
    Q_eff_bondwise(i,nearest_neighbours_index{i}) = Q_O*A*exp(-0.25*(d(i,nearest_neighbours_index{i})).^2);
    Q_eff(i)=sum(Q_eff_bondwise(i,nearest_neighbours_index{i}));
    t(nOrbitals*(i-1)+1:nOrbitals*(i),nOrbitals*(i-1)+1:nOrbitals*(i))=-eff_pot*(Q_Cu + Q_eff(i));
    
%     d_x(i,further_neighbours_index{i}) = min([abs(r(i,1)-r(further_neighbours_index{i},1)),...
%         abs(r(i,1)-r(further_neighbours_index{i},1)-N_lin),...
%         abs(r(i,1)-r(further_neighbours_index{i},1)+N_lin)],[],2);
%     d_y(i,further_neighbours_index{i}) = min([abs(r(i,2)-r(further_neighbours_index{i},2)),...
%         abs(r(i,2)-r(further_neighbours_index{i},2)-N_lin),...
%         abs(r(i,2)-r(further_neighbours_index{i},2)+N_lin)],[],2);

    
    [min_abs_d_x,min_abs_d_x_index]=min([abs(r(i,1)-r(further_neighbours_index{i},1)),...
        abs(r(i,1)-r(further_neighbours_index{i},1)-N_lin),...
        abs(r(i,1)-r(further_neighbours_index{i},1)+N_lin)],[],2);
    
    dummy_x_dist_colmns=-[(r(i,1)-r(further_neighbours_index{i},1)),...
        (r(i,1)-r(further_neighbours_index{i},1)-N_lin),...
        (r(i,1)-r(further_neighbours_index{i},1)+N_lin)];
    
    for it_min=1:length(min_abs_d_x_index)
        min_signed_d_x(it_min)=dummy_x_dist_colmns(it_min,min_abs_d_x_index(it_min));
    end
    
    for it_neigh_index=1:numel(further_neighbours_index{i})
        d_x(i,further_neighbours_index{i}(it_neigh_index)) = min_signed_d_x(it_neigh_index);
    end
    
    [min_abs_d_y,min_abs_d_y_index]=min([abs(r(i,2)-r(further_neighbours_index{i},2)),...
        abs(r(i,2)-r(further_neighbours_index{i},2)-N_lin),...
        abs(r(i,2)-r(further_neighbours_index{i},2)+N_lin)],[],2);
    
    dummy_y_dist_colmns=-[(r(i,2)-r(further_neighbours_index{i},2)),...
        (r(i,2)-r(further_neighbours_index{i},2)-N_lin),...
        (r(i,2)-r(further_neighbours_index{i},2)+N_lin)];
    
    for it_min=1:length(min_abs_d_y_index)
        min_signed_d_y(it_min)=dummy_y_dist_colmns(it_min,min_abs_d_y_index(it_min));
    end
    
    for it_neigh_index=1:numel(further_neighbours_index{i})
        d_y(i,further_neighbours_index{i}(it_neigh_index)) = min_signed_d_y(it_neigh_index);
    end
    
    d(i,further_neighbours_index{i}) = (d_x(i,further_neighbours_index{i}).^2 + d_y(i,further_neighbours_index{i}).^2).^(0.5);
    
    
    
    if bilayer_int
        [min_abs_d_x,min_abs_d_x_index]=min([abs(r(i,1)-r(intlayerNN_index{i},1)),...
            abs(r(i,1)-r(intlayerNN_index{i},1)-N_lin),...
            abs(r(i,1)-r(intlayerNN_index{i},1)+N_lin)],[],2);
        
        dummy_x_dist_colmns=-[(r(i,1)-r(intlayerNN_index{i},1)),...
            (r(i,1)-r(intlayerNN_index{i},1)-N_lin),...
            (r(i,1)-r(intlayerNN_index{i},1)+N_lin)];
        
        for it_min=1:length(min_abs_d_x_index)
            min_signed_d_x(it_min)=dummy_x_dist_colmns(it_min,min_abs_d_x_index(it_min));
        end
        
        for it_neigh_index=1:numel(intlayerNN_index{i})
            d_x(i,intlayerNN_index{i}(it_neigh_index)) = min_signed_d_x(it_neigh_index);
        end
        
        [min_abs_d_y,min_abs_d_y_index]=min([abs(r(i,2)-r(intlayerNN_index{i},2)),...
            abs(r(i,2)-r(intlayerNN_index{i},2)-N_lin),...
            abs(r(i,2)-r(intlayerNN_index{i},2)+N_lin)],[],2);
        
        dummy_y_dist_colmns=-[(r(i,2)-r(intlayerNN_index{i},2)),...
            (r(i,2)-r(intlayerNN_index{i},2)-N_lin),...
            (r(i,2)-r(intlayerNN_index{i},2)+N_lin)];
        
        for it_min=1:length(min_abs_d_y_index)
            min_signed_d_y(it_min)=dummy_y_dist_colmns(it_min,min_abs_d_y_index(it_min));
        end
        
        for it_neigh_index=1:numel(intlayerNN_index{i})
            d_y(i,intlayerNN_index{i}(it_neigh_index)) = min_signed_d_y(it_neigh_index);
        end
        
        d(i,intlayerNN_index{i}) = (d_x(i,intlayerNN_index{i}).^2 + d_y(i,intlayerNN_index{i}).^2).^(0.5);
    end
    
    if ~isempty(t_r_ref)
        for it_orb=1:nOrbitals
            t(nOrbitals*(i-1)+it_orb,nOrbitals*(further_neighbours_index{i}-1)+it_orb) = griddata(t_r_ref(:,1),t_r_ref(:,2),(10^(-3))*t_r_ref(:,3),...
                d_x(i,further_neighbours_index{i}),d_y(i,further_neighbours_index{i}),'natural');
%             if bilayer_int
%                 if it_orb==1
%                     t(nOrbitals*(i-1)+it_orb,nOrbitals*(further_neighbours_index{i}-1)+it_orb+1) = 0.08*t(nOrbitals*(i-1)+it_orb,nOrbitals*(further_neighbours_index{i}-1)+it_orb);
%                 end
%                 if it_orb==2
%                     t(nOrbitals*(i-1)+it_orb,nOrbitals*(further_neighbours_index{i}-1)+it_orb-1) = 0.08*t(nOrbitals*(i-1)+it_orb,nOrbitals*(further_neighbours_index{i}-1)+it_orb);
%                 end
%             end
        end
        if figg            
            hold on;
            for j=1:size(further_neighbours{i},1)
                tempcol=abs(t(nOrbitals*(i-1)+1,nOrbitals*(further_neighbours_index{i}(j)-1)+1));
                
                tempx_dir=abs(r(i,1)-further_neighbours{i}(j,1))';
                tempx_righ=abs(r(i,1)-further_neighbours{i}(j,1)-N_lin)';
                tempx_lef=abs(r(i,1)-further_neighbours{i}(j,1)+N_lin)';
                                
                tempy_dir=abs(r(i,2)-further_neighbours{i}(j,2))';
                tempy_righ=abs(r(i,2)-further_neighbours{i}(j,2)-N_lin)';
                tempy_lef=abs(r(i,2)-further_neighbours{i}(j,2)+N_lin)';
                
                tempx=(min([tempx_dir;...
                    tempx_righ;...
                    tempx_lef]));                
                tempy=(min([tempy_dir;...
                    tempy_righ;...
                    tempy_lef]));
                
                if and(tempx==tempx_dir,...
                        tempy==tempy_dir)
                    
                    line([r(i,1),r(i,1)+d_x(i,further_neighbours_index{i}(j))],...
                        [r(i,2),r(i,2)+d_y(i,further_neighbours_index{i}(j))],...
                        'Color',[0.90-2*tempcol,0.90-2*tempcol,0.90-6*tempcol],'Linewidth',40*tempcol);
                    hold on;
                else
                    line([r(i,1),r(i,1)+d_x(i,further_neighbours_index{i}(j))],...
                        [r(i,2),r(i,2)+d_y(i,further_neighbours_index{i}(j))],...
                        'Color',[0.90-2*tempcol,0.90-2*tempcol,0.90-6*tempcol],'LineStyle',':','Linewidth',40*tempcol);
                    hold on;
%                 elseif and(tempx==tempx_dir,...
%                         tempy==tempy_righ)
%                     line([r(i,1),r(i,1)+(further_neighbours{i}(j,1)-r(i,1))],...
%                         [r(i,2),r(i,2)+(further_neighbours{i}(j,2)-r(i,2))],...
%                         'Color',[0.90-2*tempcol,0.90-4*tempcol,0.90-6*tempcol],'Linewidth',30*tempcol,'LineStyle','-');
%                     hold on;
%                 elseif and(tempx==tempx_dir,...
%                         tempy==tempy_lef)
%                     line([r(i,1),r(i,1)+(further_neighbours{i}(j,1)-r(i,1))],...
%                         [r(i,2),r(i,2)+(further_neighbours{i}(j,2)-r(i,2))],...
%                         'Color',[0.90-2*tempcol,0.90-4*tempcol,0.90-6*tempcol],'Linewidth',30*tempcol,'LineStyle','-');
%                     hold on;
%                 elseif and(tempx==tempx_righ,...
%                         tempy==tempy_dir)
%                     line([r(i,1),r(i,1)+(further_neighbours{i}(j,1)-r(i,1))],...
%                         [r(i,2),r(i,2)+(further_neighbours{i}(j,2)-r(i,2))],...
%                         'Color',[0.90-2*tempcol,0.90-4*tempcol,0.90-6*tempcol],'Linewidth',30*tempcol,'LineStyle','-');
%                     hold on;
%                 elseif and(tempx==tempx_righ,...
%                         tempy==tempy_righ)
%                     line([r(i,1),r(i,1)+(further_neighbours{i}(j,1)-r(i,1))],...
%                         [r(i,2),r(i,2)+(further_neighbours{i}(j,2)-r(i,2))],...
%                         'Color',[0.90-2*tempcol,0.90-4*tempcol,0.90-6*tempcol],'Linewidth',30*tempcol,'LineStyle','-');
%                     hold on;
%                 elseif and(tempx==tempx_righ,...
%                         tempy==tempy_lef)
%                     line([r(i,1),r(i,1)+(further_neighbours{i}(j,1)-r(i,1))],...
%                         [r(i,2),r(i,2)+(further_neighbours{i}(j,2)-r(i,2))],...
%                         'Color',[0.90-2*tempcol,0.90-4*tempcol,0.90-6*tempcol],'Linewidth',30*tempcol,'LineStyle','-');
%                     hold on;
%                 elseif and(tempx==tempx_lef,...
%                         tempy==tempy_dir)
%                     line([r(i,1),r(i,1)+(further_neighbours{i}(j,1)-r(i,1))],...
%                         [r(i,2),r(i,2)+(further_neighbours{i}(j,2)-r(i,2))],...
%                         'Color',[0.90-2*tempcol,0.90-4*tempcol,0.90-6*tempcol],'Linewidth',30*tempcol,'LineStyle','-');
%                     hold on;
%                 elseif and(tempx==tempx_lef,...
%                         tempy==tempy_righ)
%                     line([r(i,1),r(i,1)+(further_neighbours{i}(j,1)-r(i,1))],...
%                         [r(i,2),r(i,2)+(further_neighbours{i}(j,2)-r(i,2))],...
%                         'Color',[0.90-2*tempcol,0.90-4*tempcol,0.90-6*tempcol],'Linewidth',30*tempcol,'LineStyle','-');
%                     hold on;
%                 elseif and(tempx==tempx_lef,...
%                         tempy==tempy_lef)
%                     line([r(i,1),r(i,1)+(further_neighbours{i}(j,1)-r(i,1))],...
%                         [r(i,2),r(i,2)+(further_neighbours{i}(j,2)-r(i,2))],...
%                         'Color',[0.90-2*tempcol,0.90-4*tempcol,0.90-6*tempcol],'Linewidth',30*tempcol,'LineStyle','-');
%                     hold on;                    
                end
            end
        end
        if bilayer_int
           
%             if it_orb==1
                t(nOrbitals*(i-1)+1,nOrbitals*(intlayerNN_index{i}-1)+2)=griddata(perpTBlatticeVector(:,1),perpTBlatticeVector(:,2),perpTB,d_x(i,intlayerNN_index{i}),d_y(i,intlayerNN_index{i}),'natural');
%             end;
%             if it_orb==2
                t(nOrbitals*(i-1)+2,nOrbitals*(intlayerNN_index{i}-1)+1)=t(nOrbitals*(i-1)+1,nOrbitals*(intlayerNN_index{i}-1)+2);
%             end;
                t(nOrbitals*(i-1)+1,nOrbitals*(i-1)+2)=+0.08*0.15*(1); %%%%%%% This takes into account the same (x,y) interlayer hopping
                t(nOrbitals*(i-1)+2,nOrbitals*(i-1)+1)=+0.08*0.15*(1);
        end

    else
        for it_orb=1:nOrbitals
            t(nOrbitals*(i-1)+it_orb,nOrbitals*(further_neighbours_index{i}-1)+it_orb) = -1/exp(-1)*exp(-(d(i,further_neighbours_index{i})).^2);
        end
    end
    
    
    
    if int_soc
        for it_orb=1:nOrbitals
            for it_neighb_ind=1:numel(nearest_neighbours_index{i})
                if and( atan2(d_y(i,nearest_neighbours_index{i}(it_neighb_ind)),d_x(i,nearest_neighbours_index{i}(it_neighb_ind))) > -3*pi/4,...
                        atan2(d_y(i,nearest_neighbours_index{i}(it_neighb_ind)),d_x(i,nearest_neighbours_index{i}(it_neighb_ind))) < -pi/4 )
                    
                    H_soc1(nOrbitals*(i-1)+it_orb,nOrbitals*(nearest_neighbours_index{i}(it_neighb_ind)-1)+it_orb) = 0.03*0.15*1i*(-1)^(it_orb-1);
                    
                elseif and( atan2(d_y(i,nearest_neighbours_index{i}(it_neighb_ind)),d_x(i,nearest_neighbours_index{i}(it_neighb_ind))) > pi/4,...
                        atan2(d_y(i,nearest_neighbours_index{i}(it_neighb_ind)),d_x(i,nearest_neighbours_index{i}(it_neighb_ind))) < 3*pi/4 )
                    
                    H_soc1(nOrbitals*(i-1)+it_orb,nOrbitals*(nearest_neighbours_index{i}(it_neighb_ind)-1)+it_orb) = 0.03*0.15*(-1i)*(-1)^(it_orb-1);
                    
                elseif and( atan2(d_y(i,nearest_neighbours_index{i}(it_neighb_ind)),d_x(i,nearest_neighbours_index{i}(it_neighb_ind))) > -pi/4,...
                        atan2(d_y(i,nearest_neighbours_index{i}(it_neighb_ind)),d_x(i,nearest_neighbours_index{i}(it_neighb_ind))) < pi/4 )
                    
                    H_soc2(nOrbitals*(i-1)+it_orb,nOrbitals*(nearest_neighbours_index{i}(it_neighb_ind)-1)+it_orb) = 0.03*0.15*(-1i)*(-1)^(it_orb-1); %flipping sign (1i to -1i), SO sign check, 30nov2020
                    
                elseif or( atan2(d_y(i,nearest_neighbours_index{i}(it_neighb_ind)),d_x(i,nearest_neighbours_index{i}(it_neighb_ind))) > 3*pi/4,...
                        atan2(d_y(i,nearest_neighbours_index{i}(it_neighb_ind)),d_x(i,nearest_neighbours_index{i}(it_neighb_ind))) < -3*pi/4 ) 
                    
                    H_soc2(nOrbitals*(i-1)+it_orb,nOrbitals*(nearest_neighbours_index{i}(it_neighb_ind)-1)+it_orb) = 0.03*0.15*(1i)*(-1)^(it_orb-1); %flipping sign (-1i to 1i), SO sign check, 30nov2020
               
                end
            end
        end
    end
        
        
end

% if bilayer_int
%         for i = 1: size(r,1)
%             t(nOrbitals*(i-1)+1,nOrbitals*(i-1)+2)=+0.08*0.15*(1);
%             t(nOrbitals*(i-1)+2,nOrbitals*(i-1)+1)=+0.08*0.15*(1);
%         end
% end

% signed_d_x=zeros(size(r,1));
% signed_d_y=zeros(size(r,1));
% for i=1:size(r,1)
%     for j=1:nearest_neighbours_index{i})
%         if ((r(i,1)-r(j,1))^2 + (r(i,1)-r(j,1))^2)^0.5 < NN_cut_off,   
%        
%     end
% end
    
    

if figg
    hold on;
    t_diag_dummy=diag(t);
    scatter(r(:,1),r(:,2),4*4*25+(t_diag_dummy(1:nOrbitals:end)/10).^2,-t_diag_dummy(1:nOrbitals:end)/abs(eff_pot),'filled');
    colorbar;
    c=colorbar;
    c.Label.String='Q';
    blue_red_map(gcf);ca=caxis;mx_ca=max(abs(ca));caxis([-mx_ca,mx_ca]);
    
    if further_N_cut_off < 2
        title('');
%         title({' ',[repmat(' ',1,40),'NN cut off ',num2str(NN_cut_off),', NNN cut off ',num2str(further_N_cut_off),', effective onsite charge',' (N=',num2str(N_lin),', disloc. len.=',num2str(N_lin^2-size(r,1)),')'],' '});
    else
        title('');
%         title({' ',[repmat(' ',1,40),'NN cut off ',num2str(NN_cut_off),', NNNN cut off ',num2str(further_N_cut_off),', effective onsite charge',' (N=',num2str(N_lin),', disloc. len.=',num2str(N_lin^2-size(r,1)),')'],' '});
        
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    dummygca=gca;dummygca.FontSize=25;
    axis equal;
%     axis off;
    xlim([-N_lin/2,N_lin/2]);
    ylim([-N_lin/2,N_lin/2]);
    box on;
    
    if exist('chg_map_filepath','var')
        if further_N_cut_off < 2
            print_pdf([chg_map_filepath,'/','NN cut off ',num2str(NN_cut_off),',NNN cut off ',num2str(further_N_cut_off),',effective onsite charge',' ( N=',num2str(N_lin),'disloc len=',num2str(N_lin^2-size(r,1)),' )','.pdf']);
        else
            print_pdf([chg_map_filepath,'/','NN cut off ',num2str(NN_cut_off),',NNNN cut off ',num2str(further_N_cut_off),',effective onsite charge',' ( N=',num2str(N_lin),'disloc len=',num2str(N_lin^2-size(r,1)),' )','.pdf']);
        end
    end
end
end