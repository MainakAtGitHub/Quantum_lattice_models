function [ rowColIdx ] = findRowOrColumnInMatrix( m, v )
%FINDROWORCOLUMNINMATRIX given a matrix and a vector returns the rows or
%column where this vector can be found
% the function finds a vector as either a row or a column of the supplied
% matrix if is contains the vector. The vector have to fill up entire row
% or column.
% $Author: T. Giselsson $

if nargin<2
    % create test
    m=ones(10,3);
    v=[3 4 5];
    m(5,:)=v;
    m(9,:)=v;
    m=m';
    m([3 4 5])=v;
    m=m';
end

msize=size(m);
assert(length(msize)==2,'bad matrix dimension');
v=v(:)';
vsize=size(v,2);
assert(msize(1)==vsize || msize(2)==vsize,'matrix and vector needs to share one dimension');

if msize(2)==vsize
    m=m';
end
msize=size(m);

vectorIdx=strfind(reshape(m,1,[]),v);
[i,j]=ind2sub(msize,vectorIdx);

rowColIdx=[];
for k=1:length(i)
    if i(k)~=1
        continue;
    end
    rowColIdx=[rowColIdx j(k)];
end

end