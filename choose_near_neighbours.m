function [near_neighbours,near_neighbour_index] = choose_near_neighbours(r,N_lin,d_cut_off,plot_near_neighbours,inv_con_len)
%%%%%%%temporary for reduced xbox calculation
N_linx=N_lin;%29.606;%N_lin;%
N_liny=N_lin;
%%%%%%%temporary for reduced xbox calculation
if nargin < 2
    N_lin=0;
end

if nargin < 3
    d_cut_off=1.01;
end

if nargin < 4
    plot_near_neighbours=false;
end

if nargin < 5
    inv_con_len=2;
end

near_neighbours{size(r,1),1}=[];
near_neighbour_index{size(r,1),1}=[];

for i = 1:size(r,1)    
     dist_all_neighbours = sqrt((min([abs(r(:,1)-r(i,1))';abs(r(:,1)-r(i,1)-N_linx)';...
         abs(r(:,1)-r(i,1)+N_linx)'])).^2+(min([abs(r(:,2)-r(i,2))';...
         abs(r(:,2)-r(i,2)-N_lin)';abs(r(:,2)-r(i,2)+N_lin)'])).^2);
     near_neighbour_index{i} = find(dist_all_neighbours<d_cut_off & dist_all_neighbours~=0)';
     near_neighbours{i}=[r(near_neighbour_index{i},1),r(near_neighbour_index{i},2)];
end



if plot_near_neighbours
    figure;
    for i = 1:size(r,1)
        for j=1:size(near_neighbours{i},1)
            if and((min([abs(r(i,1)-near_neighbours{i}(j,1))';...
                    abs(r(i,1)-near_neighbours{i}(j,1)-N_lin)';...
                    abs(r(i,1)-near_neighbours{i}(j,1)+N_lin)']))==abs(r(i,1)-near_neighbours{i}(j,1))',...
                    (min([abs(r(i,2)-near_neighbours{i}(j,2))';...
                    abs(r(i,2)-near_neighbours{i}(j,2)-N_lin)';...
                    abs(r(i,2)-near_neighbours{i}(j,2)+N_lin)']))==abs(r(i,2)-near_neighbours{i}(j,2))')
               plot([r(i,1),(near_neighbours{i}(j,1)+inv_con_len*r(i,1))/(inv_con_len+1)],...
                   [r(i,2),(near_neighbours{i}(j,2)+inv_con_len*r(i,2))/(inv_con_len+1)],...
                    'Color',[0, 0, 1],'Linewidth',1.5);
%                     'Color',[0.75, 0.75, 0],'Linewidth',1.5);
            else
                inv_con_len_indirect=1.5*N_lin;
                plot([r(i,1),(near_neighbours{i}(j,1)+inv_con_len_indirect*r(i,1))/(inv_con_len_indirect+1)],...
                    [r(i,2),(near_neighbours{i}(j,2)+inv_con_len_indirect*r(i,2))/(inv_con_len_indirect+1)],...
                    'Color',[1, 0, 0],'Linewidth',1);
%                     'Color',[0.8500, 0.3250, 0.0980],'Linewidth',1);
            end
            hold on;
        end
    end
    title('Neighbours');
    axis equal;
    set(gcf, 'Position', get(0, 'Screensize'));
end
