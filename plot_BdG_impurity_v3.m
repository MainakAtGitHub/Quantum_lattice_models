function p=plot_BdG_impurity_v3(inputfile)
load(inputfile,'-mat');
fig1=figure;
subplot(2,2,1); plot(nAcc); title('nAcc');% axis('square');
subplot(2,2,2); plot(muAcc); title('mu');% axis('square');
subplot(2,2,3); plot(abs(deltaMaxAcc)); title('max(deltaNN)'); %axis('square');
subplot(2,2,4); plot(abs(deltaMinAcc)); title('max(deltaNNN)'); %axis('square');

fig2=figure;% plot(deltaDiffAcc); title('Norm deltaDiff'); axis('square');
% Create semilogy
axes1 = axes('Parent',fig2,'YScale','log','YMinorTick','on',...
    'YMinorGrid','on');
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
semilogy(deltaDiffAcc);

% Create title
title('Norm deltaDiff');

% Create xlabel
xlabel('n');

% Create ylabel
ylabel('|\Delta_{n+1}-\Delta_n|');

