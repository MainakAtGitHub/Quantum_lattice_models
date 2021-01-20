function t = lattice_translation_v2(N, TBparameters, latticeVectors)
% Forms t_ij from t_i0 where t_ij represents hopping/interaction from cell
% 0 to cell i. output is N^2*nOrbitals X N^2*nOrbitals matrix. Matrix
% elements are arranged on the 2D lattice such  that columns are scanned
% first. eg. (1,1)->(1,2)-> .... (1,N) (2,1)->(2,2)....
% Row gives hopping from all other lattice sites to the Row lattice site.


% to be done: generalize to systems that are not NxN but N(1)xN(2), i.e.
% rectangular by default
if numel(N)==1
     N=[N N];
end
nOrbitals = size(TBparameters,1);
maxHop = max(abs(latticeVectors(:,1:2)),[],1);
minSystemSize = 2*maxHop + 1;
% no special treatment for hoppings that extend differently in x and y directions
% just make the system big enougth in both directions
if sum(N < minSystemSize) >0  
    % that is if either xsize is less than minsyssize or ysize, sum is (-Mainak)
    % working as "or" here (-Mainak)
    error(['System size must be greater than or equal',num2str(minSystemSize)]);
end
% to be done 1)flag for boundary conditions
% 2) sparse matrix algebra switch
t = zeros(N(1)*N(2)*nOrbitals,N(1)*N(2)*nOrbitals);
for ix = 1:N(1)
    for iy = 1:N(2)
        for jx = 1:N(1)
            for jy = 1:N(2)
                diffX = ix - jx;
                diffY = iy - jy;
                if (abs(diffX) <  N(1) - maxHop(1)) && (abs(diffY) <  N(2) - maxHop(2))
                    % direct hop
                    ind = find((latticeVectors(:,1) == diffX) & (latticeVectors(:,2) == diffY));
                    if ~isempty(ind)
                        iRange = (((ix-1)*N(2) + iy - 1)*nOrbitals + 1):((ix-1)*N(2) + iy)*nOrbitals;
                        jRange = (((jx-1)*N(2) + jy - 1)*nOrbitals + 1):((jx-1)*N(2) + jy)*nOrbitals;
                        t(iRange, jRange) = hoppings(iRange,jRange);%TBparameters(:,:,ind);
                    end
                elseif (abs(diffX) <  N(1) - maxHop(1)) && (abs(diffY) >=  N(2) - maxHop(2))
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX) & (latticeVectors(:,2) == diffY - sign(diffY)*N(2)));
                    if ~isempty(ind)
                        iRange = (((ix-1)*N(2) + iy - 1)*nOrbitals + 1):((ix-1)*N(2) + iy)*nOrbitals;
                        jRange = (((jx-1)*N(2) + jy - 1)*nOrbitals + 1):((jx-1)*N(2) + jy)*nOrbitals;
                        t(iRange, jRange) = hoppings(iRange,jRange);%TBparameters(:,:,ind);
                    end
                elseif (abs(diffY) <  N(2) - maxHop(2)) && (abs(diffX) >=  N(1) - maxHop(1))
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N(1)) & (latticeVectors(:,2) == diffY));
                    if ~isempty(ind)
                        iRange = (((ix-1)*N(2) + iy - 1)*nOrbitals + 1):((ix-1)*N(2) + iy)*nOrbitals;
                        jRange = (((jx-1)*N(2) + jy - 1)*nOrbitals + 1):((jx-1)*N(2) + jy)*nOrbitals;
                        t(iRange, jRange) = hoppings(iRange,jRange);%TBparameters(:,:,ind);
                    end
                else
                   ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N(1)) & (latticeVectors(:,2) == diffY - sign(diffY)*N(2))); 
                   if ~isempty(ind)
                        iRange = (((ix-1)*N(2) + iy - 1)*nOrbitals + 1):((ix-1)*N(2) + iy)*nOrbitals;
                        jRange = (((jx-1)*N(2) + jy - 1)*nOrbitals + 1):((jx-1)*N(2)+ jy)*nOrbitals;
                        t(iRange, jRange) = hoppings(iRange,jRange);%TBparameters(:,:,ind);
                    end
                end
            end
        end
    end
end