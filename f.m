function d1=f(E,a,kx,ky,e)
% assuming a rectangular grid in kx and ky with the corresponding energies E
% we calculate the contribution to
% the density of states from every triangle and add up all of them.
d1=zeros(1,length(e));
% d2=d1;
%fprintf(1,['\nof ',num2str(size(kx,1)),' are done:  ']);
for i=1:size(kx,1)
%       Printing Progress to the Command Window
       %testoutput(i,i-1);
       i1=i;
       if i1==size(kx,1)
           i1=1;
       end;
       % use a second variable for the contribution of one row of triangles
       % to reduce the error from the additions
       tmp=0;
    for j=1:size(ky,2)
        j1=j;
        if j1==size(ky,2)
            j1=1;
        end;
     E1=[E(i1,j1);E(i1,j1+1);E(i1+1,j1)]; %triangle without point i+1,j+1
     k1=[kx(i1,j1),ky(i1,j1);kx(i1,j1+1),ky(i1,j1+1);kx(i1+1,j1),ky(i1+1,j1)];
     a1=[a(i1,j1);a(i1,j1+1);a(i1+1,j1)];
     r1=triF(E1,a1,k1,e);
     E1=[E(i1+1,j1+1);E(i1,j1+1);E(i1+1,j1)]; %triangle without point i,j
     k1=[kx(i1+1,j1+1),ky(i1+1,j1+1);kx(i1,j1+1),ky(i1,j1+1);kx(i1+1,j1),ky(i1+1,j1)];
     a1=[a(i1+1,j1+1);a(i1,j1+1);a(i1+1,j1)];
     r2=triF(E1,a1,k1,e);
%      E1=[E(i,j);E(i+1,j);E(i+1,j+1)]; % triangle without point i,j+1
%      k1=[kx(i,j),ky(i,j);kx(i+1,j),ky(i+1,j);kx(i+1,j+1),ky(i+1,j+1)];
%      r3=triD(E1,k1,e);
%      E1=[E(i+1,j+1);E(i,j+1);E(i,j)]; % triangle without point i+1,j
%      k1=[kx(i+1,j+1),ky(i+1,j+1);kx(i,j+1),ky(i,j+1);kx(i,j),ky(i,j)];
%      r4=triD(E1,k1,e);
     % in order to minimize the errors from the triangle discritization we
     % calculate the contribution twice with different triangles and weight
     % them equally (apparently there is no difference, so we save the
     % computing time)
     tmp=tmp+(r1+r2);
%      d2=d2+abs(r1+r2-r3-r4);
    end
    d1=d1+tmp;
end
% d2par
d1=d1/((2*pi)^2);
