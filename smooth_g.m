function smooth=smooth_g(data,sigma,N)
if nargin < 3
    N=ceil(3*sigma);
end;
if nargin < 2
    sigma=1;
end;
%=conv2(Im,gaussian2d(N,sig),'same');
g=gaussian3d(2*N,sigma);
smoothc=convn(data,g);
szd=size(data)
smooth=smoothc(N+1:szd(1)+N,N+1:szd(2)+N,N+1:szd(3)+N);
size(smooth)
end
function f=gaussian3d(N,sigma)
if nargin <2
    sigma=1;
end;
  % N is grid size, sigma speaks for itself
 [x y z]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2), round(-N/2):round(N/2));
 f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2)-z.^2/(2*sigma^2));
 f=f./sum(f(:))
end