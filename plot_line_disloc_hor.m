function tt=plot_line_disloc_hor(NN,full_dis_length)
% NN=12;
if nargin<1
    NN=12;
end
if nargin<2
    full_dis_length=11;
end
points=-NN:NN;
dis_length=(full_dis_length-1)/2;
[xxx,yyy]=meshgrid(points,points);
x=xxx;
y=yyy;
grid_points=[x(:),y(:)];
dislocn=zeros(size(grid_points));
grid_points_dis=zeros(size(grid_points));

figure;
scatter(grid_points(:,1),grid_points(:,2),'o','b');
hold on;
i_iter=0;
for i=1:length(grid_points(:,1)) 
     i_iter =i_iter+1;
    if ~and(grid_points(i,2)==NN-floor(NN),any(-dis_length:dis_length==grid_points(i,1))) 
        if ~(grid_points(i,2)==NN-floor(NN))
            if grid_points(i,1)>0
                dislocn(i,:)=line_disloc_hor(grid_points(i,:),[dis_length,(full_dis_length/2)-floor(full_dis_length/2)],1);% + line_disloc_hor(grid_points(i,:),[-5,0],1))/2;
                grid_points_dis(i_iter,:)=grid_points(i,:)+ dislocn(i,:); 
            elseif grid_points(i,1)==0
                dislocn(i,:)=(line_disloc_hor(grid_points(i,:),[dis_length,(full_dis_length/2)-floor(full_dis_length/2)],1)+ line_disloc_hor(grid_points(i,:),[-dis_length,(full_dis_length/2)-floor(full_dis_length/2)],1))/2;
                grid_points_dis(i_iter,:)=grid_points(i,:)+ dislocn(i,:);
            else
                dislocn(i,:)=line_disloc_hor(grid_points(i,:),[-dis_length,(full_dis_length/2)-floor(full_dis_length/2)],1);% + line_disloc_hor(grid_points(i,:),[5,0],1))/2;
                grid_points_dis(i_iter,:)=grid_points(i,:)+ dislocn(i,:); 
            end
                % scatter(grid_points(i,1),grid_points(i,2),'b');
        else
            dislocn(i,:)=0;
            grid_points_dis(i_iter,:)=grid_points(i,:)+dislocn(i,:);
        end
        scatter(grid_points_dis(i_iter,1),grid_points_dis(i_iter,2),'.','r');
        hold on;
    end
    if and(grid_points(i,2)==NN-floor(NN),any(-dis_length:dis_length==grid_points(i,1)))
        i_iter = i_iter-1;
%         grid_points_dis(i_iter,:)=[];
    end
end
grid_points_dis = grid_points_dis(1:i_iter,:);
axis equal;
% grid_points_dis(288,:)=[];
save('data_line_disloc_hor.mat','grid_points','grid_points_dis');
tt=hoppings(grid_points_dis,2*NN+1);