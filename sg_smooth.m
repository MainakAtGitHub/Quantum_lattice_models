function sdata=sg_smooth(data,order)
% smoothing according to Savitzky–Golay filter 
% http://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter_for_smoothing_and_differentiation
if nargin <2
    order=1;
end;
switch order
    case 1
        % smoothing according to Savitzky–Golay filter
        % http://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter_for_smoothing_and_differentiation
        coeff=[-3 12 17 12 -3];
    case {2,3,4,5,6,7,8,9,10}
        % moving average
        % http://en.wikipedia.org/wiki/Moving_average
        n0=order;
        coeff=ones(1,n0);
    case 12
        coeff=[-2 3 6 7 6 3 -2];
    case 13
        coeff=[-21 14 39 54 59 54 39 14 -21];
    case 14
        coeff=[5 -30 75 131 75 -30 5];
    case 15
        coeff=[15 -55 30 135 179 135 30 -55 15];
end
if order==0
    % no smoothing
    sdata=data;
else
    lcoeff=floor(length(coeff)/2);
    snorm=sum(coeff);
    ldat=size(data);
    sdata=data;
    for n=lcoeff+1:ldat(2)-lcoeff
        sdata(:,n)=sum(data(:,n-lcoeff:n+lcoeff).*repmat(coeff,ldat(1),1),2)/snorm;
    end;
end
