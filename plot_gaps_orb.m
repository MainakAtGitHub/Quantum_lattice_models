function p=plot_gaps_orb(inputfile,all,band)
if nargin < 1
    inputfile='/home/UFAD/mainak.pal/Desktop/summer_2019_office/andreas/BdG/Tutorial/f21/no_imp/5band_dd_pppp_xy_0000/input_5band_SC.txt'
end
if nargin < 2
    all=false;
end;
if nargin < 3
    band=false;
end;
    % load the input file to set the variables, gave up the old .mat file
    % format
     read_input_file=inputfile;
     read_input;
     read_input_file
     
     close all
%nOrbitals = 2;
for which_orbital=1:nOrbitals
    if band && (which_orbital > 1)
        continue
    end
    for which_orbital2=1:nOrbitals
        if ~all
            if ~(which_orbital2==which_orbital)
                continue
            end
        end
        figure(which_orbital+(which_orbital2-1)*nOrbitals);
        xlabel('k_x/\pi')
        ylabel('k_y/\pi')
        if band
            title(['\Delta_{band}'])
        else
            title(['\Delta_{',num2str(which_orbital), num2str(which_orbital2),'}']);
        end
        axis equal
    end
end;

N = 100;
k_range =-pi:pi/N:pi;
[KX,KY] = meshgrid(k_range);
klist = [KX(:),KY(:)];
[eigv,eigvc]=bands_precalc_klist(inputfile,klist);
delta_max=0;
for j=1:nOrbitals
    clearvars fermi_k_path_dummy posn numb patch;
    eigv_reshaped = reshape(eigv(:,j),length(k_range),length(k_range));
    fermi_k_path_dummy=contourc(k_range,k_range,eigv_reshaped,[0,0])';
    if size(fermi_k_path_dummy,1)> 0
        %fermi_k_path_dummy = contourc(k_range,k_range,eigv_reshaped(:,:,j),[0,0])';
        posn(1) = 1;
        numb(1) = fermi_k_path_dummy(posn(1),2);
        fermi_k_path{j}= fermi_k_path_dummy(posn(1)+1:posn(1)+numb(1),:);
        i = 1;
        while true
            posn(i+1) = posn(i)+numb(i)+1;
            if posn(i+1) > length(fermi_k_path_dummy(:,1))
                break
            end
            numb(i+1) = fermi_k_path_dummy(posn(i+1),2);
            patch = fermi_k_path_dummy(posn(i+1)+1:posn(i+1)+numb(i+1),:);
            fermi_k_path{j} = [fermi_k_path{j};patch];
            i=i+1;
        end
        delta_mat = real(Deltak_klist(inputfile,fermi_k_path{j}(:,:)));
        for which_orbital=1:nOrbitals
            if band && (which_orbital > 1)
                continue
            end
                
            for which_orbital2=1:nOrbitals
                if ~all
                    if ~(which_orbital2==which_orbital)
                        continue
                    end
                end
                figure(which_orbital+(which_orbital2-1)*nOrbitals);
                hold on
                if ~band
                    delta = delta_mat(:,which_orbital,which_orbital2);
                else
                    %code for band transformation
                    [~,kSpaceEigenVectorsNormal]=bands_precalc_klist(inputfile,fermi_k_path{j}(:,:));
                    delta = zeros(1,length(fermi_k_path{j}(:,1)));
                    for ii = 1:length(fermi_k_path{j}(:,1))
                       delta_band = squeeze(kSpaceEigenVectorsNormal(ii,:,:))*squeeze(delta_mat(ii,:,:))*squeeze(kSpaceEigenVectorsNormal(ii,:,:))';
                       delta(ii) = delta_band(j,j);
                    end   
                end
               
                delta_max=max([delta_max,max(abs(delta(:)))]);
                % stem3(fermi_k_path{j}(:,1)/pi,fermi_k_path{j}(:,2)/pi,delta,'.');
                delta = real(delta);
                scatter3(fermi_k_path{j}(:,1)/pi,fermi_k_path{j}(:,2)/pi,delta,20,delta);
                %if j==nOrbitals
                caxis([-delta_max delta_max]);
                colorbar;
                %end
            end
        end
    end
end 

