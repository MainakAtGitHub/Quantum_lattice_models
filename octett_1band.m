function [o,figure1]=octett_1band(kx,ky,energy,figure1)
% visualize the octett model given the kx, ky points for the ends of the
% "banananas" as calculated using "banana_1band.m"
% definition of the onset point as in Hanaguri et al. 
qx=-ky(:)/pi;
qy=-kx(:)/pi;
% give back a struct with all the q-vectors
o.q1=[-qx-qx,qy-qy];
o.q2=[-qy-qx,qx-qy];
o.q3=[-qy-qx,-qx-qy];
o.q4=[-qx-qx,-qy-qy];
o.q5=[qx-qx,-qy-qy];
o.q6=[qy-qx,-qx-qy];
o.q7=[qy-qx,qx-qy];
% make a fancy plot
if nargin <4
    figure1=figure;
end;
if nargin <4
scatter(o.q1(:,1),o.q1(:,2),100,energy,'o','filled','DisplayName','q_1')
findaxis=findall(figure1,'type','axes');
for n=1:numel(findaxis)
    hold(findaxis(n),'on');
end;
scatter(o.q2(:,1),o.q2(:,2),100,energy,'v','filled','DisplayName','q_2')
scatter(o.q3(:,1),o.q3(:,2),100,energy,'^','filled','DisplayName','q_3')
scatter(o.q4(:,1),o.q4(:,2),100,energy,'s','filled','DisplayName','q_4')
scatter(o.q5(:,1),o.q5(:,2),100,energy,'p','filled','DisplayName','q_5')
scatter(o.q6(:,1),o.q6(:,2),100,energy,'h','filled','DisplayName','q_6')
scatter(o.q7(:,1),o.q7(:,2),100,energy,'d','filled','DisplayName','q_7')
else
    findaxis=findall(figure1,'type','axes');
for n=1:numel(findaxis)
    hold(findaxis(n),'on');
end;
    scatter3(o.q1(:,1),o.q1(:,2),10,100,'o','filled','DisplayName','q_1')
scatter3(o.q2(:,1),o.q2(:,2),10,100,'v','filled','DisplayName','q_2')
scatter3(o.q3(:,1),o.q3(:,2),10,100,'^','filled','DisplayName','q_3')
scatter3(o.q4(:,1),o.q4(:,2),10,100,'s','filled','DisplayName','q_4')
scatter3(o.q5(:,1),o.q5(:,2),10,100,'p','filled','DisplayName','q_5')
scatter3(o.q6(:,1),o.q6(:,2),10,100,'h','filled','DisplayName','q_6')
scatter3(o.q7(:,1),o.q7(:,2),10,100,'d','filled','DisplayName','q_7')
end
for n=1:numel(findaxis)
    hold(findaxis(n),'off');
end;
if nargin <4
    legend1 = legend('show');
    % Create colorbar
    colorbar;
end
xlabel('q_x/\pi');
ylabel('q_y/\pi');
% c = 1:numel(t);      %# colors
% h = surface([x(:), x(:)], [y(:), y(:)], [z(:), z(:)], ...
%     [c(:), c(:)], 'EdgeColor','flat', 'FaceColor','none');
% colormap( jet(numel(t)) )

% plot also a 2D plot as in Hanaguri talk
fig2=figure;
% plot(energy,sqrt(o.q1(:,1).^2+o.q1(:,2).^2),'o','DisplayName','q_1');
% hold on
% plot(energy,sqrt(o.q2(:,1).^2+o.q2(:,2).^2),'v','DisplayName','q_2');
% plot(energy,sqrt(o.q3(:,1).^2+o.q3(:,2).^2),'^','DisplayName','q_3');
% plot(energy,sqrt(o.q4(:,1).^2+o.q4(:,2).^2),'s','DisplayName','q_4');
% plot(energy,sqrt(o.q5(:,1).^2+o.q5(:,2).^2),'p','DisplayName','q_5');
% plot(energy,sqrt(o.q6(:,1).^2+o.q6(:,2).^2),'h','DisplayName','q_6');
% plot(energy,sqrt(o.q7(:,1).^2+o.q7(:,2).^2),'d','DisplayName','q_7');

plot(energy,sqrt(o.q1(:,1).^2+o.q1(:,2).^2),'k','DisplayName','q_1');
hold on
plot(energy,sqrt(o.q2(:,1).^2+o.q2(:,2).^2),'b','DisplayName','q_2');
plot(energy,sqrt(o.q3(:,1).^2+o.q3(:,2).^2)','r','DisplayName','q_3');
plot(energy,sqrt(o.q4(:,1).^2+o.q4(:,2).^2),'g','DisplayName','q_4');
plot(energy,sqrt(o.q5(:,1).^2+o.q5(:,2).^2),'--k','DisplayName','q_5');
plot(energy,sqrt(o.q6(:,1).^2+o.q6(:,2).^2),'--r','DisplayName','q_6');
plot(energy,sqrt(o.q7(:,1).^2+o.q7(:,2).^2),'--b','DisplayName','q_7');

 legend('show');
xlabel('E (meV)');
ylabel('|q|/\pi');
