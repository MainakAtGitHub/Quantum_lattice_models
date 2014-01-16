% inputs
nOrbitals = 10;
M = 5;
ita = .001;
E = -.0084;
xGridRange = -80:60;
yGridRange = -60:80;
zGridRange = [0 21 4 22];
load lattice_greens_supercell_FeSe_N_15_M_9_U_0955_Vimp_5_E_minPt0084.mat
load wannier_FeSe_4d_matrix_v2


nBands = size(latticeGreens,1);
N = sqrt(nBands/nOrbitals);

% Local greens function
% Wannier vector
RDiscrete = [40 40 80];
shift = [51 51 41];
sizeWannier = [101 101 81];
fileName = ['supercell_local_ldos_FeSe_U_0955_Vimp_5','_N_',num2str(N),'_M_',num2str(M),'_E_',num2str(E),'_ita_',num2str(ita)];
for zGridPoint = zGridRange
    loacalLdos = zeros(length(xGridRange),length(yGridRange));
    countLoopX = 0;
    for xGridPoint = xGridRange
        countLoopX = countLoopX + 1;
        countLoopY = 0;
        for yGridPoint = yGridRange
            countLoopY = countLoopY + 1;
            r = [xGridPoint, yGridPoint, zGridPoint];
            wAcc = []; 
            for l = -(ceil(N/2)-1):(ceil(N/2)-1)
                for m = -(ceil(N/2)-1):(ceil(N/2)-1)
                    for orbital = 1:nOrbitals
                        R = [l m 0];
                        latticeVector = RDiscrete.*R;
                        wannierArgument = r - latticeVector;                       
                        shiftedArgument = wannierArgument + shift; % translate wannier origin                                     
                        % check whether this argument is in range or not
                        if (((shiftedArgument) <= sizeWannier) & ((shiftedArgument) >= [1 1 1]))
                            % yes in range, now find the value
                            %wannierValue = wannierI(shiftedArgument);
                            w = wannierValues(shiftedArgument(1),shiftedArgument(2),shiftedArgument(3),orbital);
                        else
                            % not in range, set it to zero
                            w = 0;
                        end
                        wAcc = [wAcc; w];
                    end
                end
            end
            loacalLdos(countLoopX,countLoopY) = (-1/pi)*imag(wAcc'*(latticeGreens*wAcc));
        end
    end
    fileName = ['supercell_spatial_ldos_FeSe_N_15_M_5_U_0955_Vimp_50_E_-0.0084_ita_0.001', '_z_', num2str(zGridPoint),'.mat'];
    save([fileName2,'_z_',num2str(zGridPoint),'.mat'],'loacalLdos','xGridRange','yGridRange');
end

