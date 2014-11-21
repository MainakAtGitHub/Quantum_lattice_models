function t = lattice_single(N, Impparameters, ImpVectors, center)
% Forms t_ij from t_i0 where t_ij represents hopping/interaction from cell
% 0 to cell i. output is N^2*nOrbitals X N^2*nOrbitals matrix. Matrix
% elements are arranged on the 2D lattice such  that columns are scanned
% first. eg. (1,1)->(1,2)-> .... (1,N) (2,1)->(2,2)....
% Row gives hopping from all other lattice sites to the Row lattice site.

nOrbitals = size(Impparameters,1);
maxHop = max(max(abs(ImpVectors)));
minSystemSize = 2*maxHop + 1;
if N < minSystemSize
    error(['System size must be greater than or equal',num2str(minSystemSize)]);
end
t = zeros(N^2*nOrbitals,N^2*nOrbitals);
% run over all ImpVectors
szImpVectors=size(ImpVectors);
% center position of impurity
% make the impurity vectors only 2 dimensional
if szImpVectors(2)==6
    ImpVectors=ImpVectors(:,[1:2,4:5]);
end;
for n=1:szImpVectors(1)
    i=ImpVectors(n,1:2)+center;
    j=ImpVectors(n,3:4)+center;
    % periodic hop if needed
    i=mod(i-1,N)+1;
    j=mod(j-1,N)+1;
    iRange = (((i(1)-1)*N + i(2) - 1)*nOrbitals + 1):((i(1)-1)*N + i(2))*nOrbitals;
    jRange = (((j(1)-1)*N + j(2) - 1)*nOrbitals + 1):((j(1)-1)*N + j(2))*nOrbitals;
     t(iRange, jRange) = Impparameters(:,:,n);
end;
% 
% for ix = 1:N
%     for iy = 1:N
%         for jx = 1:N
%             for jy = 1:N
%                 diffX = ix - jx;
%                 diffY = iy - jy;
%                 if (abs(diffX) <  N - maxHop) && (abs(diffY) <  N - maxHop)
%                     % direct hop
%                     ind = find((ImpVectors(:,1) == diffX) & (ImpVectors(:,2) == diffY));
%                     if ~isempty(ind)
%                         iRange = (((ix-1)*N + iy - 1)*nOrbitals + 1):((ix-1)*N + iy)*nOrbitals;
%                         jRange = (((jx-1)*N + jy - 1)*nOrbitals + 1):((jx-1)*N + jy)*nOrbitals;
%                         t(iRange, jRange) = Impparameters(:,:,ind);
%                     end
%                 elseif (abs(diffX) <  N - maxHop) && (abs(diffY) >=  N - maxHop)
%                     % periodic hop
%                     ind = find((ImpVectors(:,1) == diffX) & (ImpVectors(:,2) == diffY - sign(diffY)*N));
%                     if ~isempty(ind)
%                         iRange = (((ix-1)*N + iy - 1)*nOrbitals + 1):((ix-1)*N + iy)*nOrbitals;
%                         jRange = (((jx-1)*N + jy - 1)*nOrbitals + 1):((jx-1)*N + jy)*nOrbitals;
%                         t(iRange, jRange) = Impparameters(:,:,ind);
%                     end
%                 elseif (abs(diffY) <  N - maxHop) && (abs(diffX) >=  N - maxHop)
%                     % periodic hop
%                     ind = find((ImpVectors(:,1) == diffX - sign(diffX)*N) & (ImpVectors(:,2) == diffY));
%                     if ~isempty(ind)
%                         iRange = (((ix-1)*N + iy - 1)*nOrbitals + 1):((ix-1)*N + iy)*nOrbitals;
%                         jRange = (((jx-1)*N + jy - 1)*nOrbitals + 1):((jx-1)*N + jy)*nOrbitals;
%                         t(iRange, jRange) = Impparameters(:,:,ind);
%                     end
%                 else
%                     % periodic hop in x and y direction
%                    ind = find((ImpVectors(:,1) == diffX - sign(diffX)*N) & (ImpVectors(:,2) == diffY - sign(diffY)*N)); 
%                    if ~isempty(ind)
%                         iRange = (((ix-1)*N + iy - 1)*nOrbitals + 1):((ix-1)*N + iy)*nOrbitals;
%                         jRange = (((jx-1)*N + jy - 1)*nOrbitals + 1):((jx-1)*N + jy)*nOrbitals;
%                         t(iRange, jRange) = Impparameters(:,:,ind);
%                     end
%                 end
%             end
%         end
%     end
% end