function smooth=smooth_g(data,sigma,periodic)
% smooth_g    gaussian smoothing
% smooth=smooth_g(data) smoothing of data (3d array) with standart deviation sigma=1
% and non-periodic borders
%
% somooth=smooth_g(data,sigma,periodic) smoothing with given standart deviation in
% sigma and periodicity given in periodic (default periodic=[ 0 0 0] , no
% periodicity in all three directions

% do smoothing over a range of 6 sigma
N=ceil(6*max(sigma));
if nargin < 2
    sigma=1;
end;
% generate smoothing matrix
g=gaussian3d(N,sigma);
% number for extension of array
M=(size(g,1)-1)/2;
% extend data periodically dimensios given in variable "periodic"
data=padarray(data,periodic*M,'circular','both');
szd=size(data);
nonperiodic=(periodic==0);
% just copy the edges to extend the data in all other directions
data=padarray(data,nonperiodic*M,'replicate');
% keep track of extension of matrix
periodic=periodic+nonperiodic;
% convolution of data with smoothing matrix
smoothc=convn(data,g);
clear data;
% remove the additional data from smoothing the edges
smooth=smoothc(M*(1+periodic(1))+1:szd(1)+M*(1+periodic(1)),M*(1+periodic(2))+1:szd(2)+M*(1+periodic(2)),M*(1+periodic(3))+1:szd(3)+M*(1+periodic(3)));
% size(smooth) % debug output
end
function f=gaussian3d(N,sigma)
if nargin <2
    sigma=1;
end;
  % N is grid size, sigma speaks for itself
 [x y z]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2), round(-N/2):round(N/2));
 % generate some gaussian matrix with standart deviation sigma
 if numel(sigma)==1
     sigma=ones(1,3)*sigma;
 end;
 f=exp(-x.^2/(2*sigma(1)^2)-y.^2/(2*sigma(2)^2)-z.^2/(2*sigma(3)^2));
 f=f./sum(f(:));
end