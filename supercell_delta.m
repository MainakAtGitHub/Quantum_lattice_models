function [deltaSuper,superDeltaVectors] = supercell_delta(nOrbitals, delta, maxHop)

%maxHop = 2;
N = sqrt(size(delta,1)/nOrbitals);
deltaSuper = zeros(N^2*nOrbitals,N^2*nOrbitals,9);
superDeltaVectors = [0 0; 1 0; -1 0; 0 1; 0 -1; 1 1; -1 1; -1 -1; 1 -1];

for ix = 1:N
    for iy = 1:N
        for jx = 1:N
            for jy = 1:N
                diffX = ix - jx;
                diffY = iy - jy;
                [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, [ix iy], [jx jy]);
                if (abs(diffX) <  N - maxHop) && (abs(diffY) <  N - maxHop)
                    deltaSuper(iRange, jRange, 1) = delta(iRange, jRange);
                elseif (diffX <=  -N + maxHop) && (abs(diffY) <  N - maxHop)
                    deltaSuper(iRange, jRange, 2) = delta(iRange, jRange);
                elseif (diffX >=  N - maxHop) && (abs(diffY) <  N - maxHop)
                    deltaSuper(iRange, jRange, 3) = delta(iRange, jRange);
                elseif (diffY <=  -N + maxHop) && (abs(diffX) <  N - maxHop)
                    deltaSuper(iRange, jRange, 4) = delta(iRange, jRange);
                elseif (diffY >=  N - maxHop) && (abs(diffX) <  N - maxHop)
                    deltaSuper(iRange, jRange, 5) = delta(iRange, jRange);
                elseif (diffX <=  -N + maxHop) && (diffY <=  -N + maxHop)
                    deltaSuper(iRange, jRange, 6) = delta(iRange, jRange);
                elseif (diffX >=  N - maxHop) && (diffY <=  -N + maxHop)
                    deltaSuper(iRange, jRange, 7) = delta(iRange, jRange);
                elseif (diffX >=  N - maxHop) && (diffY >=  N - maxHop)
                    deltaSuper(iRange, jRange, 8) = delta(iRange, jRange);
                elseif (diffX <=  -N + maxHop) && (diffY >=  N - maxHop)
                    deltaSuper(iRange, jRange, 9) = delta(iRange, jRange);
                end
            end
        end
    end
end
                    
                    
                    
                    
                    