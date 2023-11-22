function [ output_args ] = plot_triplet_order( inputfile, which_del ) % whichDelta can be 1,2,3,4 for UpUp, UpDown,..
 
        read_input_file=inputfile;
        read_input;
        read_input_file
        load(pos_file,'-mat');
        load(BdGfileName);
        [filepath,name,ext]=fileparts(BdGfileName);
        if exist('tripletOrderFile','var')
            load(tripletOrderFile); % a file like Gamma files, to multiply 1,i,-1,-i to the right, up, left, down bonds to get triplet order
            
            tripletOrderMatrix=lattice_translation(N, tripletOrder, tripletOrderVecs);
            if which_del==1
                delta=delta(1:size(delta,1)/2, 1:size(delta,1)/2);
            elseif which_del==2
                delta=delta(1:size(delta,1)/2, (size(delta,1)/2 + 1):(size(delta,1)/2 + size(delta,1)/2));
            elseif which_del==3
                delta=(delta((size(delta,1)/2 + 1):(size(delta,1)/2 + size(delta,1)/2), 1:size(delta,1)/2));
            elseif which_del==4
                delta=delta((size(delta,1)/2 + 1):(size(delta,1)/2 + size(delta,1)/2), (size(delta,1)/2 + 1):(size(delta,1)/2 + size(delta,1)/2));
            end
            tripletOrderDelta=tripletOrderMatrix.*delta;
            figure; imagesc(real(tripletOrderDelta));colorbar;axis equal;blue_red_map(gcf);
            figure; imagesc(imag(tripletOrderDelta));colorbar;axis equal;blue_red_map(gcf);
            figure; imagesc(abs(tripletOrderDelta));colorbar;axis equal;bluemap(gcf);
            max(max(abs(real(tripletOrderDelta))))
            min(min(abs(real(tripletOrderDelta))))
            mean(mean(abs(real(tripletOrderDelta))))
            max(max(abs(imag(tripletOrderDelta))))
            min(min(abs(imag(tripletOrderDelta))))
            mean(mean(abs(imag(tripletOrderDelta))))
% N = sqrt(size(delta,1)/nOrbitals);
% if exist('latticeVectorsSC','var')
%     nUnitCellsDelta = size(latticeVectorsSC,1);
% else
%     latticeVectorsSC=Gammafull.latt;
%     nUnitCellsDelta = size(Gammafull.latt,1);
% end
% 
% 
% deltamap = zeros(N, N, nOrbitals,nOrbitals);
% % read out delta and store it into map of matrices
% % to do generalize to non-rectangular systems
% for N1=1:N
%     for N2=1:N
%         jCell = [N1 N2];
%         for i = 1:nUnitCellsDelta
%             iCell = jCell + latticeVectorsSC(i,:);
%             cellvector=periodic_latticevectors(iCell,N);
%             [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, cellvector, jCell);
%             position=jCell+latticeVectorsSC(i,:)/2;
%             % now add the contribution to the 4 neighbored lattice points
%             pos(1,:)=periodic_latticevectors(floor(position),N);
%             pos(2,:)=periodic_latticevectors(ceil(position),N);
%             pos(3,:)=periodic_latticevectors([ceil(position(1)),floor(position(2))],N);
%             pos(4,:)=periodic_latticevectors([floor(position(1)),ceil(position(2))],N);
%             for ps=1:4
%                 deltamap(pos(ps,1),pos(ps,2),:,:) =deltamap(pos(ps,1),pos(ps,2),:,:)+shiftdim(delta(iRange, jRange).*conj(delta(iRange, jRange)),-2);
%             end
%         end
%     end
% end
        end
end

