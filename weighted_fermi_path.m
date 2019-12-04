function p=weighted_fermi_path(inputfile,all,band)
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
     
     N = 200;
k_range =-pi:pi/N:pi;
[KX,KY] = meshgrid(k_range);
klist = [KX(:),KY(:)];
[eigv,eigvc]=bands_precalc_klist(inputfile,klist);
delta_max=0;
figure();
% colorlist = [1,0,0; 0,0.502,0.502; 0.5020,0.5020,0; 1,0,1; 0,1,1];
colorlist = [1,0,0; 0,1,0; 0,0,1; 1,0,1; 0,1,1];
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
        [~,kSpaceEigenVectorsNormal]=bands_precalc_klist(inputfile,fermi_k_path{j}(:,:));
        weight = zeros(1,length(fermi_k_path{j}(:,1)));
        orb_color = zeros(length(fermi_k_path{j}(:,1)),3);
        orb_numb = zeros(1,length(fermi_k_path{j}(:,1)));
        for ii = 1:length(fermi_k_path{j}(:,1))
            jth_vec = kSpaceEigenVectorsNormal(ii,j,:);
            weight(ii) = max(abs(jth_vec).^2)/sum(abs(jth_vec).^2);
            orb_numb(ii) = find(abs(jth_vec) == max(abs(jth_vec)));
            orb_color(ii,:) = colorlist(orb_numb(ii),:);
        end        
%         orb1=orb_color==1;
%         w1=weight(orb1);
%         kx1=fermi_k_path{j}(orb1,1)/pi;
%         ky1=fermi_k_path{j}(orb1,2)/pi;
%                 
%         scatter3(kx1,ky1,w1,'MarkerEdgeColor',[0 0.75 0.75],'MarkerFaceColor',[0 .75 .75]);

        scatter3(fermi_k_path{j}(:,1)/pi,fermi_k_path{j}(:,2)/pi,weight,20,orb_color,'filled');
%         legend(orb_color(ii,:));
        colorbar;
        hold on
    end
end

