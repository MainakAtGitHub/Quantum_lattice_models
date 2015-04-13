function [t, superLatticeVectors] = supercell_hoppings(N, TBparameters, latticeVectors)

nOrbitals = size(TBparameters,1);
maxHop = max(max(abs(latticeVectors)));
t = zeros(N^2*nOrbitals,N^2*nOrbitals,9);
superLatticeVectors = [0 0; 1 0; -1 0; 0 1; 0 -1; 1 1; -1 1; -1 -1; 1 -1];

for ix = 1:N
    for iy = 1:N
        for jx = 1:N
            for jy = 1:N
                diffX = ix - jx;
                diffY = iy - jy;
                [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, [ix iy], [jx jy]);
                if (abs(diffX) <  N - maxHop) && (abs(diffY) <  N - maxHop)
                    % direct hop
                    ind = find((latticeVectors(:,1) == diffX) & (latticeVectors(:,2) == diffY));
                    if ~isempty(ind)
                        t(iRange, jRange,1) = TBparameters(:,:,ind);
                    end
                elseif (diffX <=  -N + maxHop) && (abs(diffY) <  N - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N) & (latticeVectors(:,2) == diffY));
                    if ~isempty(ind)
                        t(iRange, jRange,2) = TBparameters(:,:,ind);
                    end
                elseif (diffX >=  N - maxHop) && (abs(diffY) <  N - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N) & (latticeVectors(:,2) == diffY));
                    if ~isempty(ind)
                        t(iRange, jRange,3) = TBparameters(:,:,ind);
                    end
                elseif (diffY <=  -N + maxHop) && (abs(diffX) <  N - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX) & (latticeVectors(:,2) == diffY - sign(diffY)*N));
                    if ~isempty(ind)
                        t(iRange, jRange,4) = TBparameters(:,:,ind);
                    end
                elseif (diffY >=  N - maxHop) && (abs(diffX) <  N - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX) & (latticeVectors(:,2) == diffY - sign(diffY)*N));
                    if ~isempty(ind)
                        t(iRange, jRange,5) = TBparameters(:,:,ind);
                    end
                elseif (diffX <=  -N + maxHop) && (diffY <=  -N + maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N) & (latticeVectors(:,2) == diffY - sign(diffY)*N));
                    if ~isempty(ind)
                        t(iRange, jRange,6) = TBparameters(:,:,ind);
                    end
                elseif (diffX >=  N - maxHop) && (diffY <=  -N + maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N) & (latticeVectors(:,2) == diffY - sign(diffY)*N));
                    if ~isempty(ind)
                        t(iRange, jRange,7) = TBparameters(:,:,ind);
                    end
                elseif (diffX >=  N - maxHop) && (diffY >=  N - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N) & (latticeVectors(:,2) == diffY - sign(diffY)*N));
                    if ~isempty(ind)
                        t(iRange, jRange,8) = TBparameters(:,:,ind);
                    end
                elseif (diffX <=  -N + maxHop) && (diffY >=  N - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N) & (latticeVectors(:,2) == diffY - sign(diffY)*N));
                    if ~isempty(ind)
                        t(iRange, jRange,9) = TBparameters(:,:,ind);
                    end
                end
            end
        end
    end
end
                    
                    
