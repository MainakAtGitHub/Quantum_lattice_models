function [o,figure1]=plot_octett_1band(kx,ky,energy,figure1,symm,h)
if nargin <5
    symm=false;
end;
if nargin < 6
    h=0.00002;
end;
% visualize the octett model given the kx, ky points for the ends of the
% "banananas" as calculated using "banana_1band.m"
% definition of the onset point as in Hanaguri et al.
[o,energy]=octett_1band(kx,ky,energy,symm);
% rotate the vectors according to the symmetry
fillstring='filled'
%fillstring='empty'
% make a fancy plot
newfig=true;
if nargin <4
    figure1=figure;
elseif isempty(figure1)
        figure1=figure;
else
    figure1=figure(figure1);
    newfig=false;
end;
shapes='ov^sphd';
col=[.6 .6 1];
col=[1 0 0];
lwth=2;
if newfig
scatter(o.q1(:,1),o.q1(:,2),h,energy,shapes(1),'filled','DisplayName','q_1','LineWidth',lwth)
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
sizedot=100;
absolut=sign(energy);
scatter3(o.q1(:,1),o.q1(:,2),h*ones(length(energy),1),sizedot,h*ones(length(energy),1),shapes(1),'DisplayName','q_1','LineWidth',lwth,'MarkerEdgeColor',col)
scatter3(o.q2(:,1),o.q2(:,2),h*ones(length(energy),1),sizedot,h*ones(length(energy),1),shapes(1),'DisplayName','q_2','LineWidth',lwth,'MarkerEdgeColor',col)
scatter3(o.q3(:,1),o.q3(:,2),h*ones(length(energy),1),sizedot,h*ones(length(energy),1),shapes(1),'DisplayName','q_3','LineWidth',lwth,'MarkerEdgeColor',col)
scatter3(o.q4(:,1),o.q4(:,2),h*ones(length(energy),1),sizedot,h*ones(length(energy),1),shapes(1),'DisplayName','q_4','LineWidth',lwth,'MarkerEdgeColor',col)
scatter3(o.q5(:,1),o.q5(:,2),h*ones(length(energy),1),sizedot,h*ones(length(energy),1),shapes(1),'DisplayName','q_5','LineWidth',lwth,'MarkerEdgeColor',col)
scatter3(o.q6(:,1),o.q6(:,2),h*ones(length(energy),1),sizedot,h*ones(length(energy),1),shapes(1),'DisplayName','q_6','LineWidth',lwth,'MarkerEdgeColor',col)
scatter3(o.q7(:,1),o.q7(:,2),h*ones(length(energy),1),sizedot,h*ones(length(energy),1),shapes(1),'DisplayName','q_7','LineWidth',lwth,'MarkerEdgeColor',col)
end
for n=1:numel(findaxis)
    hold(findaxis(n),'off');
end;
axis square
xlim([-2 2]);
ylim([-2 2]);
% annotation(figure1,'textbox',...
%     [0.6 0.1 0.3 0.08],...
%     'String',{['E= ',num2str(energy(1)*1000),' meV']},...
%     'FontSize',14,...
%     'FitBoxToText','off', 'EdgeColor','none','TextColor','red');

if newfig
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
if ~symm
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
end;