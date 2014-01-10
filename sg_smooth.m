function sdata=sg_smooth(data,order)
% smoothing according to Savitzky–Golay filter default: window size 7
% (quadratic)
% http://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter_for_smoothing_and_differentiation
if nargin <2
    order=3;
end;
if order < 0
    % moving average
    % http://en.wikipedia.org/wiki/Moving_average
    coeff=ones(1,-order*2+1);
end
if order > 0
    % smoothing according to Savitzky–Golay filter
    % http://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter_for_smoothing_and_differentiation
    coeff=conv_coeff(2*order+1);
end
if order==0
    % no smoothing
    sdata=data;
else
    ldat=size(data);
	% remove spikes
    for da=1:ldat(1)
        data(da,:)=remove_spike(data(da,:));
    end;
    lcoeff=floor(length(coeff)/2);
    snorm=sum(coeff);
    sdata=data;
    for n=lcoeff+1:ldat(2)-lcoeff
        sdata(:,n)=sum(data(:,n-lcoeff:n+lcoeff).*repmat(coeff,ldat(1),1),2)/snorm;
    end;
end
function c=conv_coeff(n)
for a=1:n
    c(a)=(3*n^2-7-20*(a-ceil(n/2))^2)/4/(n*(n^2-4))*3;
end;
function data=remove_spike(data)
pos=find(abs(data)==max(abs(data)));
ldata=length(data);
if (pos(1)>2) && (pos(1)<ldata-1)
    % ignore spikes on the edges
if (abs(data(pos))>5*abs(data(pos(1)+1)) && abs(data(pos(1)))>5*abs(data(pos(1)-1)))
	data(pos)=0.5*(data(pos+1)+data(pos-1));
	disp(['Data point ',num2str(pos),'  corrected']);
	data=remove_spike(data);
end;
end;