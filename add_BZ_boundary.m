function fulldata=add_BZ_boundary(data)
% assuming a data in a rectangular matrix (in a BZ), this function adds
% data from the boundary such that it is suitable for plotting and using in
% numerical integration (via singular_quad)

fulldata=[data,data(:,1);...
          data(1,:),data(1,1)];
% dummy test
%fulldata=data;