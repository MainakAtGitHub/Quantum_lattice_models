function [tcdelta]=fit_gap_tc(gap_set,T_set, tcmax)
% from some data points [T,Delta] fit the BCS result of the gap equation
% to extract T and Delta
if nargin <3
    tcmax=max(T_set)/2;
end;
load('mean_field.csv');
options = optimset;
options.TolX=1e-8;
options.TolFun=1e-7;
[tcdelta,fval,exitflag,output] = fminsearch(@(TcDelta)calc_weight(TcDelta,gap_set,T_set,mean_field),[tcmax/2,max(gap_set)],options);

end
function weight=calc_weight(TcDelta,gap_set,T_set,mean_field)
gap=bcs_gap(T_set/TcDelta(1),mean_field)*TcDelta(2);
weight=sum((gap-gap_set).^2)+(TcDelta(1)>max(gap_set(:)))*(TcDelta(1)-max(gap_set(:)))^2+(TcDelta(1)<0)+(TcDelta(2)<0);
end
function gap=bcs_gap(t,mean_field)
small=t<1;
gap=0*t;
% here do an interpolation
gap(small) = interp1(mean_field(:,1),mean_field(:,2),t(small));
%gap(small)=tanh(1.82*(1.018*(1-t(small)).^(0.51)));
end