function [deltaSuper,superDeltaVectors] = supercell_delta(nOrbitals, delta, maxHop,spin_and_nambu,N)

if nargin < 4
    spin_and_nambu=0;
end
%maxHop = 2;
% try
%     read_input_file=inputfile;
%     read_input;
%     read_input_file
%     if numel(N)==1
%         N=[N,N];
%     end
% catch
%     N = sqrt(size(delta,1)/nOrbitals);
%     N =[N,N];
% end

% if ~exist('spin_and_nambu','var')
%     spin_and_nambu=0;
% end
if nargin < 5
    N = sqrt(size(delta,1)/nOrbitals);
end
    if numel(N)==1
        N=[N,N];
    end

if ~spin_and_nambu
    deltaSuper = zeros(N(1)*N(2)*nOrbitals,N(1)*N(2)*nOrbitals,9);
else
    deltaSuper1 = zeros(N(1)*N(2)*nOrbitals,N(1)*N(2)*nOrbitals,9);
    deltaSuper2 = zeros(N(1)*N(2)*nOrbitals,N(1)*N(2)*nOrbitals,9);
    deltaSuper3 = zeros(N(1)*N(2)*nOrbitals,N(1)*N(2)*nOrbitals,9);
    deltaSuper4 = zeros(N(1)*N(2)*nOrbitals,N(1)*N(2)*nOrbitals,9);
end
superDeltaVectors = [0 0; 1 0; -1 0; 0 1; 0 -1; 1 1; -1 1; -1 -1; 1 -1];

if spin_and_nambu
        delta1= delta(1:size(delta,1)/2,1:size(delta,1)/2);
        delta2= delta(1:size(delta,1)/2,size(delta,1)/2+1:size(delta,1));
        delta3= delta(size(delta,1)/2+1:size(delta,1),1:size(delta,1)/2);
        delta4= delta(size(delta,1)/2+1:size(delta,1),size(delta,1)/2+1:size(delta,1));    
end

for ix = 1:N(1)
    for iy = 1:N(2)
        for jx = 1:N(1)
            for jy = 1:N(2)
                diffX = ix - jx;
                diffY = iy - jy;
                [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, [ix iy], [jx jy]);
                if (abs(diffX) <  N(1) - maxHop) && (abs(diffY) <  N(2) - maxHop)
                    if ~spin_and_nambu
                        deltaSuper(iRange, jRange, 1) = delta(iRange, jRange);
                    else
                        deltaSuper1(iRange, jRange, 1) = delta1(iRange, jRange);
                        deltaSuper2(iRange, jRange, 1) = delta2(iRange, jRange);
                        deltaSuper3(iRange, jRange, 1) = delta3(iRange, jRange);
                        deltaSuper4(iRange, jRange, 1) = delta4(iRange, jRange);
                    end
                elseif (diffX <=  -N(1) + maxHop) && (abs(diffY) <  N(2) - maxHop)
                    if ~spin_and_nambu
                        deltaSuper(iRange, jRange, 2) = delta(iRange, jRange);
                    else
                        deltaSuper1(iRange, jRange, 2) = delta1(iRange, jRange);
                        deltaSuper2(iRange, jRange, 2) = delta2(iRange, jRange);
                        deltaSuper3(iRange, jRange, 2) = delta3(iRange, jRange);
                        deltaSuper4(iRange, jRange, 2) = delta4(iRange, jRange);
                    end
                elseif (diffX >=  N(1) - maxHop) && (abs(diffY) <  N(2) - maxHop)
                    if ~spin_and_nambu
                        deltaSuper(iRange, jRange, 3) = delta(iRange, jRange);
                    else
                        deltaSuper1(iRange, jRange, 3) = delta1(iRange, jRange);
                        deltaSuper2(iRange, jRange, 3) = delta2(iRange, jRange);
                        deltaSuper3(iRange, jRange, 3) = delta3(iRange, jRange);
                        deltaSuper4(iRange, jRange, 3) = delta4(iRange, jRange);
                    end
                elseif (diffY <=  -N(2) + maxHop) && (abs(diffX) <  N(1) - maxHop)
                    if ~spin_and_nambu
                        deltaSuper(iRange, jRange, 4) = delta(iRange, jRange);
                    else
                        deltaSuper1(iRange, jRange, 4) = delta1(iRange, jRange);
                        deltaSuper2(iRange, jRange, 4) = delta2(iRange, jRange);
                        deltaSuper3(iRange, jRange, 4) = delta3(iRange, jRange);
                        deltaSuper4(iRange, jRange, 4) = delta4(iRange, jRange);
                    end
                elseif (diffY >=  N(2) - maxHop) && (abs(diffX) <  N(1) - maxHop)
                    if ~spin_and_nambu
                        deltaSuper(iRange, jRange, 5) = delta(iRange, jRange);
                    else
                        deltaSuper1(iRange, jRange, 5) = delta1(iRange, jRange);
                        deltaSuper2(iRange, jRange, 5) = delta2(iRange, jRange);
                        deltaSuper3(iRange, jRange, 5) = delta3(iRange, jRange);
                        deltaSuper4(iRange, jRange, 5) = delta4(iRange, jRange);
                    end
                elseif (diffX <=  -N(1) + maxHop) && (diffY <=  -N(2) + maxHop)
                    if ~spin_and_nambu
                        deltaSuper(iRange, jRange, 6) = delta(iRange, jRange);
                    else
                        deltaSuper1(iRange, jRange, 6) = delta1(iRange, jRange);
                        deltaSuper2(iRange, jRange, 6) = delta2(iRange, jRange);
                        deltaSuper3(iRange, jRange, 6) = delta3(iRange, jRange);
                        deltaSuper4(iRange, jRange, 6) = delta4(iRange, jRange);
                    end
                elseif (diffX >=  N(1) - maxHop) && (diffY <=  -N(2) + maxHop)
                    if ~spin_and_nambu
                        deltaSuper(iRange, jRange, 7) = delta(iRange, jRange);
                    else
                        deltaSuper1(iRange, jRange, 7) = delta1(iRange, jRange);
                        deltaSuper2(iRange, jRange, 7) = delta2(iRange, jRange);
                        deltaSuper3(iRange, jRange, 7) = delta3(iRange, jRange);
                        deltaSuper4(iRange, jRange, 7) = delta4(iRange, jRange);
                    end
                elseif (diffX >=  N(1) - maxHop) && (diffY >=  N(2) - maxHop)
                    if ~spin_and_nambu
                        deltaSuper(iRange, jRange, 8) = delta(iRange, jRange);
                    else
                        deltaSuper1(iRange, jRange, 8) = delta1(iRange, jRange);
                        deltaSuper2(iRange, jRange, 8) = delta2(iRange, jRange);
                        deltaSuper3(iRange, jRange, 8) = delta3(iRange, jRange);
                        deltaSuper4(iRange, jRange, 8) = delta4(iRange, jRange);
                    end
                elseif (diffX <=  -N(1) + maxHop) && (diffY >=  N(2) - maxHop)
                    if ~spin_and_nambu
                        deltaSuper(iRange, jRange, 9) = delta(iRange, jRange);
                    else
                        deltaSuper1(iRange, jRange, 9) = delta1(iRange, jRange);
                        deltaSuper2(iRange, jRange, 9) = delta2(iRange, jRange);
                        deltaSuper3(iRange, jRange, 9) = delta3(iRange, jRange);
                        deltaSuper4(iRange, jRange, 9) = delta4(iRange, jRange);
                    end
                end
            end
        end
    end
end
if spin_and_nambu
    deltaSuper =[deltaSuper1,deltaSuper2;deltaSuper3,deltaSuper4];
end
                    
                    
                    
                    
                    