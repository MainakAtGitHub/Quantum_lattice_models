function t = lattice_translation(N, TBparameters, latticeVectors)
% Forms t_ij from t_i0 where t_ij represents hopping/interaction from cell
% 0 to cell i. output is N^2*nOrbitals X N^2*nOrbitals matrix. Matrix
% elements are arranged on the 2D lattice such  that columns are scanned
% first. eg. (1,1)->(1,2)-> .... (1,N) (2,1)->(2,2)....
% Row gives hopping from all other lattice sites to the Row lattice site.


% to be done: generalize to systems that are not NxN but N(1)xN(2), i.e.
% rectangular

nOrbitals = size(TBparameters,1);
maxHop = max(max(abs(latticeVectors)));
minSystemSize = 2*maxHop + 1;
if N < minSystemSize
    error(['System size must be greater than or equal',num2str(minSystemSize)]);
end
t = zeros(N^2*nOrbitals,N^2*nOrbitals);
for ix = 1:N
    for iy = 1:N
        for jx = 1:N
            for jy = 1:N
                diffX = ix - jx;
                diffY = iy - jy;
                if (abs(diffX) <  N - maxHop) && (abs(diffY) <  N - maxHop)
                    % direct hop
                    ind = find((latticeVectors(:,1) == diffX) & (latticeVectors(:,2) == diffY));
                    if ~isempty(ind)
                        iRange = (((ix-1)*N + iy - 1)*nOrbitals + 1):((ix-1)*N + iy)*nOrbitals;
                        jRange = (((jx-1)*N + jy - 1)*nOrbitals + 1):((jx-1)*N + jy)*nOrbitals;
                        t(iRange, jRange) = TBparameters(:,:,ind);
                    end
                elseif (abs(diffX) <  N - maxHop) && (abs(diffY) >=  N - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX) & (latticeVectors(:,2) == diffY - sign(diffY)*N));
                    if ~isempty(ind)
                        iRange = (((ix-1)*N + iy - 1)*nOrbitals + 1):((ix-1)*N + iy)*nOrbitals;
                        jRange = (((jx-1)*N + jy - 1)*nOrbitals + 1):((jx-1)*N + jy)*nOrbitals;
                        t(iRange, jRange) = TBparameters(:,:,ind);
                    end
                elseif (abs(diffY) <  N - maxHop) && (abs(diffX) >=  N - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N) & (latticeVectors(:,2) == diffY));
                    if ~isempty(ind)
                        iRange = (((ix-1)*N + iy - 1)*nOrbitals + 1):((ix-1)*N + iy)*nOrbitals;
                        jRange = (((jx-1)*N + jy - 1)*nOrbitals + 1):((jx-1)*N + jy)*nOrbitals;
                        t(iRange, jRange) = TBparameters(:,:,ind);
                    end
                else
                   ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N) & (latticeVectors(:,2) == diffY - sign(diffY)*N)); 
                   if ~isempty(ind)
                        iRange = (((ix-1)*N + iy - 1)*nOrbitals + 1):((ix-1)*N + iy)*nOrbitals;
                        jRange = (((jx-1)*N + jy - 1)*nOrbitals + 1):((jx-1)*N + jy)*nOrbitals;
                        t(iRange, jRange) = TBparameters(:,:,ind);
                    end
                end
            end
        end
    end
end