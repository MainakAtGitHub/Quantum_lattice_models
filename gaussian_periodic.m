function g=gaussian_periodic(inputmatrix)
h = fspecial('gaussian', 7, 1);
X=repmat(inputmatrix,3,3);
Y = filter2(h,X);
g=Y(size(inputmatrix,1)+1:2*size(inputmatrix,1),size(inputmatrix,2)+1:2*size(inputmatrix,2));
end