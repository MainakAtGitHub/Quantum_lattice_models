function [t, superLatticeVectors] = supercell_hoppings(N, TBparameters, latticeVectors)

if numel(N)==1
    N=[N,N];
end

nOrbitals = size(TBparameters,1);
maxHop = max(max(abs(latticeVectors)));
t = zeros(N(1)*N(2)*nOrbitals,N(1)*N(2)*nOrbitals,9);
superLatticeVectors = [0 0; 1 0; -1 0; 0 1; 0 -1; 1 1; -1 1; -1 -1; 1 -1];

for ix = 1:N(1)
    for iy = 1:N(2)
        for jx = 1:N(1)
            for jy = 1:N(2)
                diffX = ix - jx;
                diffY = iy - jy;
                [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, [ix iy], [jx jy]);
                if (abs(diffX) <  N(1) - maxHop) && (abs(diffY) <  N(2) - maxHop)
                    % direct hop
                    ind = find((latticeVectors(:,1) == diffX) & (latticeVectors(:,2) == diffY));
                    if ~isempty(ind)
                        t(iRange, jRange,1) = TBparameters(:,:,ind);
                    end
                elseif (diffX <=  -N(1) + maxHop) && (abs(diffY) <  N(2) - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N(1)) & (latticeVectors(:,2) == diffY));
                    if ~isempty(ind)
                        t(iRange, jRange,2) = TBparameters(:,:,ind);
                    end
                elseif (diffX >=  N(1) - maxHop) && (abs(diffY) <  N(2) - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N(1)) & (latticeVectors(:,2) == diffY));
                    if ~isempty(ind)
                        t(iRange, jRange,3) = TBparameters(:,:,ind);
                    end
                elseif (diffY <=  -N(2) + maxHop) && (abs(diffX) <  N(1) - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX) & (latticeVectors(:,2) == diffY - sign(diffY)*N(2)));
                    if ~isempty(ind)
                        t(iRange, jRange,4) = TBparameters(:,:,ind);
                    end
                elseif (diffY >=  N(2) - maxHop) && (abs(diffX) <  N(1) - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX) & (latticeVectors(:,2) == diffY - sign(diffY)*N(2)));
                    if ~isempty(ind)
                        t(iRange, jRange,5) = TBparameters(:,:,ind);
                    end
                elseif (diffX <=  -N(1) + maxHop) && (diffY <=  -N(2) + maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N(1)) & (latticeVectors(:,2) == diffY - sign(diffY)*N(2)));
                    if ~isempty(ind)
                        t(iRange, jRange,6) = TBparameters(:,:,ind);
                    end
                elseif (diffX >=  N(1) - maxHop) && (diffY <=  -N(2) + maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N(1)) & (latticeVectors(:,2) == diffY - sign(diffY)*N(2)));
                    if ~isempty(ind)
                        t(iRange, jRange,7) = TBparameters(:,:,ind);
                    end
                elseif (diffX >=  N(1) - maxHop) && (diffY >=  N(2) - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N(1)) & (latticeVectors(:,2) == diffY - sign(diffY)*N(2)));
                    if ~isempty(ind)
                        t(iRange, jRange,8) = TBparameters(:,:,ind);
                    end
                elseif (diffX <=  -N(1) + maxHop) && (diffY >=  N(2) - maxHop)
                    % periodic hop
                    ind = find((latticeVectors(:,1) == diffX - sign(diffX)*N(1)) & (latticeVectors(:,2) == diffY - sign(diffY)*N(2)));
                    if ~isempty(ind)
                        t(iRange, jRange,9) = TBparameters(:,:,ind);
                    end
                end
            end
        end
    end
end
                    
                    
