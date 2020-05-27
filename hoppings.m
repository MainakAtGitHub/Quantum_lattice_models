function t = hoppings(r,N_lin,figg) % r is a :by2 matrix of coordinates
if nargin < 2
    N_lin=0;
end
d = zeros(size(r,1));
t = zeros(size(r,1));
Q_eff_dummy=zeros(size(r,1));
Q_eff=zeros(size(r,1),1);
Q_Cu=+4;
Q_O=-2;
lmbda=1;
A=0.5*exp(1/4);
if nargin<3
    figg=false;
end
% if figg
%     figure;
% end
% figure;
% for i = 1:size(r,1)
% %     for j = 1:size(r,1)
%         [B,IX] = sort(sqrt((r(:,1)-r(i,1)).^2+(r(:,2)-r(i,2)).^2));
% %       plot([r(i,1);r(IX(1:5),1)],[r(i,2);r(IX(1:5),2)],'b');
% %         plot([r(i,1),r(IX(1),1)],[r(i,2),r(IX(1),2)]);
% %         hold on;
%         plot([r(i,1),r(IX(2),1)],[r(i,2),r(IX(2),2)],'b');
%         hold on;
%         plot([r(i,1),r(IX(3),1)],[r(i,2),r(IX(3),2)],'b');
%         hold on;
%         plot([r(i,1),r(IX(4),1)],[r(i,2),r(IX(4),2)],'b');
%         hold on;
%         plot([r(i,1),r(IX(5),1)],[r(i,2),r(IX(5),2)],'b');
%         hold on;
% end
 for i = 1:size(r,1)    
     dist_all_neighbours = sqrt((min([abs(r(:,1)-r(i,1))';abs(r(:,1)-r(i,1)-N_lin)';abs(r(:,1)-r(i,1)+N_lin)'])).^2+(min([abs(r(:,2)-r(i,2))';abs(r(:,2)-r(i,2)-N_lin)';abs(r(:,2)-r(i,2)+N_lin)'])).^2);
     near_neighbour_index = find(dist_all_neighbours<1.3 & dist_all_neighbours~=0);
     dist_near_neighbours = dist_all_neighbours(near_neighbour_index);
     for j = 1:length(near_neighbour_index)
         
         %%%*#?#******######????? FIX PERIODIC HOPPING t BELOW
         d(i,near_neighbour_index(j))= sqrt((min([abs(r(near_neighbour_index(j),1)-r(i,1)),abs(r(near_neighbour_index(j),1)-r(i,1)-N_lin),abs(r(near_neighbour_index(j),1)-r(i,1)+N_lin)]))^2+(min([abs(r(near_neighbour_index(j),2)-r(i,2)),abs(r(near_neighbour_index(j),2)-r(i,2)-N_lin),abs(r(near_neighbour_index(j),2)-r(i,2)+N_lin)]))^2); 
         t(i,near_neighbour_index(j)) = -1/exp(-1)*exp(-(d(i,near_neighbour_index(j)))^2);
%          t(i,near_neighbour_index(j)) = -1/exp(-1)*exp(-sqrt((r(near_neighbour_index(j),1)-r(i,1)).^2+(r(near_neighbour_index(j),2)-r(i,2)).^2));
%          t(i,near_neighbour_index(j)) = -1;
         Q_eff_dummy(i,near_neighbour_index(j))= Q_O*A*exp(-0.25*(d(i,near_neighbour_index(j)))^2);
         if figg
             if or(abs(r(i,1)-r(near_neighbour_index(j),1))>1,abs(r(i,2)-r(near_neighbour_index(j),2))>1)
                 pl=line([r(i,1),r(near_neighbour_index(j),1)],[r(i,2),r(near_neighbour_index(j),2)]);
%                   pl.LineWidth=2*abs(t(i,near_neighbour_index(j)));
%                  pl.Color='red';
%                  pl.LineStyle='--';
             else
                 pl=line([r(i,1),r(near_neighbour_index(j),1)],[r(i,2),r(near_neighbour_index(j),2)]);
%                   pl.LineWidth=2*abs(t(i,near_neighbour_index(j)));
%                  pl.Color='blue';
%                  pl.LineStyle=':';
             end
             hold on;
         end
         Q_eff(i)=Q_eff(i)+Q_eff_dummy(i,near_neighbour_index(j));
     end
     t(i,i)=-1*(Q_Cu + Q_eff(i));
 end
 if figg
     scatter(r(:,1),r(:,2),20+5*(Q_Cu+Q_eff),Q_Cu+Q_eff,'filled');
     for ii=1:size(r,1)
%          text(r(ii,1)+0.1,r(ii,2)+0.2,num2str(ii));
     end 
 axis equal;
 colorbar;
 end
        
        