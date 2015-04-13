function h=plot_k_space_gap(inputfile)
    % load the input file to set the variables, gave up the old .mat file
    % format
     read_input_file=inputfile;
     read_input;
     read_input_file
if (~exist('Vimp','var'))
    if abs(Vimp)>0
        disp('Warning: finite impurity potential, not homogeneous case.')
    end;
end;
load(TB_file,'-mat');
nOrbitals = size(TBparameters,1);
    load(Gamma_file,'-mat');
    load(BdGfileName,'-mat');
N = sqrt(size(delta,1)/nOrbitals);
if exist('latticeVectorsSC','var')
    nUnitCellsDelta = size(latticeVectorsSC,1);
else
    latticeVectorsSC=Gammafull.latt;
    nUnitCellsDelta = size(Gammafull.latt,1);
end
deltaCenter = zeros(nOrbitals, nOrbitals, size(latticeVectorsSC,1));
jCell = [ceil(N/2) ceil(N/2)];
for i = 1:nUnitCellsDelta
    iCell = jCell + latticeVectorsSC(i,:);
    [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, iCell, jCell);
    deltaCenter(:,:,i) = delta(iRange, jRange);
end
delta = deltaCenter;
%delta(:,:,iUnitCellDelta)
mxreal=max(abs(latticeVectorsSC(:)));
M=2*mxreal+1;
for mu=1:nOrb
    for nu=1:nOrb
realgap=zeros(M,M);
        % only works for square lattice
        for iUnitCellDelta = 1:nUnitCellsDelta
            iLatticeVectorDelta = latticeVectorsSC(iUnitCellDelta,:);
            realgap(iLatticeVectorDelta(1)+mxreal+1,iLatticeVectorDelta(2)+mxreal+1)=delta(mu,nu,iUnitCellDelta);
        end
    allgaps((mu-1)*M+1:(mu)*M,(nu-1)*M+1:(nu)*M)=real(realgap);
    end;
end;



% if nOrb>1
% xz=rem(ordering,5);
% ordering=floor(ordering/5);
% yz=rem(ordering,5);
% ordering=floor(ordering/5);
% x2y2=mod(ordering,5);
% ordering=floor(ordering/5);
% xy=rem(ordering,5);
% tickx{1+xz}= '$d_{xz}$';
%         tickx{1+yz}= '$d_{yz}$';
%         tickx{1+xy}= '$d_{xy}$';
%         tickx{1+x2y2}= '$d_{x^2-y^2}$';
%         tickx{11-xz-yz-xy-x2y2}= '$d_{z^2}$';
% lab=[0:nkIntegration:(nOrb-1)*nkIntegration]+nkIntegration/2+0.5;
% % plot all orbital contributions in one figure
    q=figure;
 %   axes1 = axes('Parent',q,'XTickLabel',tickx,'YTickLabel',tickx,...
 %       'XTick',lab,'YTick',lab,...
 %       'DataAspectRatio',[1 1 .01]);
    %view(axes1,[0.5 90]);
    nkIntegration=M;
    ylim([.5 nkIntegration*nOrb+0.5])
    xlim([.5 nkIntegration*nOrb+0.5])
%    grid(axes1,'on');
 %   hold(axes1,'all');
    % Create surf
    %surf(real(gaps),'Parent',axes1,'LineStyle','none');
    maxabsekkn=max(max((real(allgaps))));
    minabsekkn=min(min((real(allgaps))));
    global colorred
colorred=false;
mtix=10^(ceil(log(maxabsekkn)/log(10)));
if (ceil(log(maxabsekkn)/log(10))-log(maxabsekkn)/log(10)>0.7)
    ticks0=[-1:.05:1]; 
else  
    ticks0=[-1:.1:1]; 
end
ticks=ticks0*mtix;
labels = num2str(repmat(sign(ticks).*(abs(ticks)), 1, 1)', 2);
sqrtscale=true;
    if sqrtscale
    image(sign(real(allgaps)).*sqrt(abs(real(allgaps)/maxabsekkn))*128+128);
    ticks=sign(ticks0).*sqrt(abs(ticks));
    maxabsekkn=sqrt(maxabsekkn);
else
    image(sign(real(allgaps)).*(abs(real(allgaps)/maxabsekkn))*128+128);
end;
fsz=10;
copydir='';
    if ~isempty(copydir)
        pth=copydir;
    else
        [pth,BdGfileName1,~]=fileparts(BdGfileName);
    end;
set(0,'DefaultAxesFontSize',fsz)
   % imagesc(real(allgaps),'Parent',axes1);

    %max2=max(maxabsekkn,abs(minabsekkn));
    %maxabsekkn=max2;
    %minabsekkn=-max2;
    color1=[1 0 0]; % red
    %color2=[0 0 0]; % black
    color2=[1 1 1]; % white
    color3=[0 0 1]; % blue
    %if ~colorred
    input=[-10 -1 0 1 10]+.5;
    colormatrix=[color1; color1;color2;color3; color3];
    %else
    %    input=[0 1 10];
    %    colormatrix=[color2;color1; color1];
    %end;
    %if ~isempty(figure1)
    % draw a colorbar
    %     if ~colorred
    %     x=-maxabsekkn:maxabsekkn/255:maxabsekkn;
    %     else
    %if abs(maxabsekkn) > abs(minabsekkn)
    x=(maxabsekkn+minabsekkn)/2:maxabsekkn/255:maxabsekkn;
    %else
    %    x=(-maxabsekkn-minabsekkn)/2:abs(minabsekkn)/255:abs(minabsekkn);
    %    x=-x(end:-1:1);
    %end
    %    end
    r=interp1(input*maxabsekkn,colormatrix(:,1),x);
    g=interp1(input*maxabsekkn,colormatrix(:,2),x);
    b=interp1(input*maxabsekkn,colormatrix(:,3),x);
    set(q,'Colormap', [r',g',b']);
    
   % colorbar
    cb=colorbar_rwb(q,maxabsekkn,ticks,labels);
 %   label_boxes(nOrb,nkIntegration,tickx);
    string=[pth,'/',BdGfileName1,'real_space_all'];
    axis square
    print('-dpng',[string,'.png']);
%end;



nUnitCells = size(latticeVector,1);
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) = ...
TBparameters(:,:,(latticeVector(:,1)==0) & (latticeVector(:,2)==0)) - mu*eye(nOrbitals);

M=N*M;
M=100;
kx = (2*pi/M)*(0:(M - 1)) + pi/M;
ky = kx;
delKx = kx(2)-kx(1);
delKy = delKx;
%energy = linspace(firstEnergy, lastEnergy, nEnergyPoints);
kSpaceEigenValues = zeros(M, M, 2*nOrbitals);
kSpaceEigenVectors = zeros(M, M, 2*nOrbitals, 2*nOrbitals);
kSpaceEigenValuesNormal = zeros(M, M, nOrbitals);
kSpacegapall = zeros(M, M, nOrbitals,nOrbitals);
kSpaceEigenVectorsNormal = zeros(M, M, nOrbitals, nOrbitals);

% some double code from homogeneous_dos_v2.m !
for iKx = 1:M
        for iKy = 1:M
            k = [kx(iKx) ky(iKy)];
            % diagonalizing for normal state DOS
            %kSpaceHopping = zeros(nOrbitals,nOrbitals);
            %for iUnitCell = 1:nUnitCells
            %    iLatticeVector = latticeVector(iUnitCell,1:2);
            %    kSpaceHopping = kSpaceHopping + TBparameters(:,:,iUnitCell)*exp(1i*(iLatticeVector*k'));
            %end                        
          %  [eigVectorNormal eigValueNormal] = eig(kSpaceHopping);
          %  [eigValueKNormal sortingIndexNormal] = sort(real(diag(eigValueNormal)));
           % eigVectorKNormal = (eigVectorNormal(:,sortingIndexNormal))';
           % kSpaceEigenValuesNormal(iKx, iKy, :) = eigValueKNormal;
           % kSpaceEigenVectorsNormal(iKx, iKy, :,:) =  eigVectorKNormal;
            % diagonalizing for SC state DOS
            kSpaceGap =zeros(nOrbitals,nOrbitals);
            for iUnitCellDelta = 1:nUnitCellsDelta
                iLatticeVectorDelta = latticeVectorsSC(iUnitCellDelta,:);
                kSpaceGap = kSpaceGap + delta(:,:,iUnitCellDelta)*exp(1i*(iLatticeVectorDelta*k'));
            end
            kSpacegapall(iKx,iKy,:,:)=kSpaceGap(:,:);
           % kSpaceHamiltonian = [kSpaceHopping -kSpaceGap; -kSpaceGap' -kSpaceHopping];
           % [eigVector eigValue] = eig(kSpaceHamiltonian);
           % [eigValueK sortingIndex] = sort(real(diag(eigValue)));
           % eigVectorK = (eigVector(:,sortingIndex))';
           % kSpaceEigenValues(iKx, iKy, :) = eigValueK;
           % kSpaceEigenVectors(iKx, iKy, :,:) =  eigVectorK;
        end
        if  mod(iKx,10)==0   
            disp(['Done ',num2str(iKx), ' of ',num2str(M),' kx values.']);
        end;
end
for mu=1:nOrb
    for nu=1:nOrb
        % only works for square lattice
    allgaps((mu-1)*M+1:(mu)*M,(nu-1)*M+1:(nu)*M)=  fftshift(real(kSpacegapall(:,:,mu,nu)));
    end;
end;

% if nOrb>1
% xz=rem(ordering,5);
% ordering=floor(ordering/5);
% yz=rem(ordering,5);
% ordering=floor(ordering/5);
% x2y2=mod(ordering,5);
% ordering=floor(ordering/5);
% xy=rem(ordering,5);
% tickx{1+xz}= '$d_{xz}$';
%         tickx{1+yz}= '$d_{yz}$';
%         tickx{1+xy}= '$d_{xy}$';
%         tickx{1+x2y2}= '$d_{x^2-y^2}$';
%         tickx{11-xz-yz-xy-x2y2}= '$d_{z^2}$';
% lab=[0:nkIntegration:(nOrb-1)*nkIntegration]+nkIntegration/2+0.5;
% % plot all orbital contributions in one figure
    q=figure;
 %   axes1 = axes('Parent',q,'XTickLabel',tickx,'YTickLabel',tickx,...
 %       'XTick',lab,'YTick',lab,...
 %       'DataAspectRatio',[1 1 .01]);
    %view(axes1,[0.5 90]);
    nkIntegration=M;
    ylim([.5 nkIntegration*nOrb+0.5])
    xlim([.5 nkIntegration*nOrb+0.5])
%    grid(axes1,'on');
 %   hold(axes1,'all');
    % Create surf
    %surf(real(gaps),'Parent',axes1,'LineStyle','none');
    maxabsekkn=max(max((real(allgaps))));
    minabsekkn=min(min((real(allgaps))));
    global colorred
colorred=false;
mtix=10^(ceil(log(maxabsekkn)/log(10)));
if (ceil(log(maxabsekkn)/log(10))-log(maxabsekkn)/log(10)>0.7)
    ticks0=[-1:.05:1]; 
else  
    ticks0=[-1:.1:1]; 
end
ticks=ticks0*mtix;
labels = num2str(repmat(sign(ticks).*(abs(ticks)), 1, 1)', 2);
sqrtscale=true;
    if sqrtscale
    image(sign(real(allgaps)).*sqrt(abs(real(allgaps)/maxabsekkn))*128+128);
    ticks=sign(ticks0).*sqrt(abs(ticks));
    maxabsekkn=sqrt(maxabsekkn);
else
    image(sign(real(allgaps)).*(abs(real(allgaps)/maxabsekkn))*128+128);
end;
fsz=20;
copydir='';
    if ~isempty(copydir)
        pth=copydir;
    else
        [pth,BdGfileName1,~]=fileparts(BdGfileName);
    end;
set(0,'DefaultAxesFontSize',fsz)
   % imagesc(real(allgaps),'Parent',axes1);

    %max2=max(maxabsekkn,abs(minabsekkn));
    %maxabsekkn=max2;
    %minabsekkn=-max2;
    color1=[1 0 0]; % red
    %color2=[0 0 0]; % black
    color2=[1 1 1]; % white
    color3=[0 0 1]; % blue
    %if ~colorred
    input=[-10 -1 0 1 10]+.5;
    colormatrix=[color1; color1;color2;color3; color3];
    %else
    %    input=[0 1 10];
    %    colormatrix=[color2;color1; color1];
    %end;
    %if ~isempty(figure1)
    % draw a colorbar
    %     if ~colorred
    %     x=-maxabsekkn:maxabsekkn/255:maxabsekkn;
    %     else
    %if abs(maxabsekkn) > abs(minabsekkn)
    x=(maxabsekkn+minabsekkn)/2:maxabsekkn/255:maxabsekkn;
    %else
    %    x=(-maxabsekkn-minabsekkn)/2:abs(minabsekkn)/255:abs(minabsekkn);
    %    x=-x(end:-1:1);
    %end
    %    end
    r=interp1(input*maxabsekkn,colormatrix(:,1),x);
    g=interp1(input*maxabsekkn,colormatrix(:,2),x);
    b=interp1(input*maxabsekkn,colormatrix(:,3),x);
    set(q,'Colormap', [r',g',b']);
    
   % colorbar
    cb=colorbar_rwb(q,maxabsekkn,ticks,labels);
 %   label_boxes(nOrb,nkIntegration,tickx);
    string=[pth,'/',BdGfileName1,'k_space_all'];
    axis square
    print('-dpng',[string,'.png']);
%end;
