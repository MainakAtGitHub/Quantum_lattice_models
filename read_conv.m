function m=read_conv(T)
if nargin >0
N=[11 13 15];
%T='_06';
for nN=1:length(N)
    conv(:,nN)=plot_real_space(['./calc/SC_U/FeSe_N_',num2str(N(nN)),'_U_0955',T,'.mat'],'');
end;

else
    N=[7:2:25];
    for nN=1:length(N)
    conv(:,nN)=plot_real_space(['./calc/SC_U/FeSe_N_',num2str(N(nN)),'_U_0955.mat'],'');
end;
T='02';
end;

figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
%% Uncomment the following line to preserve the X-limits of the axes
 ylim(axes1,[0 max(conv(:))]);
%% Uncomment the following line to preserve the Y-limits of the axes
 xlim(axes1,[0 max(1./N.^2)]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to plot
plot1 = plot(1./N.^2,conv','Parent',axes1,'MarkerSize',15);
set(plot1(1),'Marker','x','DisplayName','on site');
set(plot1(2),'Marker','*','DisplayName','NN');
set(plot1(3),'Marker','.','DisplayName','NNN');

% Create xlabel
xlabel('1/N^2');

% Create ylabel
ylabel('\Delta [meV]');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.711704702627939 0.225909537856441 0.159232365145228 0.174655850540806]);
print_pdf(['/tmp/FeSe_conv_T',T,'.pdf']);
% give back an estimate of the gap for N-> infinity
 m=sum(repmat(N.^2.,3,1).*conv,2)/sum(N.^2)