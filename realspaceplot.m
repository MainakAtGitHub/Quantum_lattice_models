function figure1=realspaceplot(data2plot,numl,tickx,flnm,scale)
if nargin < 4
    scale='s';
end;
figure1=figure;
global colorred
colorred=false; % set colorscale
global fsz
fsz = 25; % font
 datarealmax = max(abs(real(data2plot(:)))); % for setting colormap
% if scale == 's'
%     image(sign(real(delta2Plot)).*sqrt(abs(real(delta2Plot)/deltaRealMax))*128+128);
% else
%     image(sign(real(delta2Plot)).*(abs(real(delta2Plot)/deltaRealMax))*128+128);
% end
% colorbar_rwb(figure1,deltaRealMax,ticks);
% label_boxes(nOrbitals/2,numl,tickx);

% set the scale
lm=log(datarealmax)/log(10);
mtix=10^(ceil(lm));
% do some refinement to avoid only single labels
if (ceil(lm)-lm > 0.5)
    tx=[-.5:.05:.5];
elseif    (ceil(lm)-lm > 0.2)
    tx=[-1:0.1:1];
else
    tx=[-1:0.2:1];
end;
ticks=mtix*tx; 
labels = num2str(repmat(sign(ticks).*(abs(ticks)), 1, 1)', 2);
if scale=='s'
    image(sign(real(data2plot)).*sqrt(abs(real(data2plot)/datarealmax))*128+128);
    ticks=sign(tx).*sqrt(mtix*abs(tx));
    datarealmax=sqrt(datarealmax);
else
    image(sign(real(data2plot)).*(abs(real(data2plot)/datarealmax))*128+128);
end;
colorbar_rwb(figure1,datarealmax,ticks,labels);
label_boxes(numel(tickx),numl,tickx);
if nargin >4
if isunix
    k = strfind(flnm, '/');
    if isempty(k)
        stringp=['/tmp/',flnm(1:length(flnm)-4),'.pdf']
    else
        stringp=['/tmp/',flnm(k(length(k))+1:length(flnm)-4),'.pdf']
    end;
    %stringp=['/tmp/Gapfunction_realspace_imag.pdf'];
    print_pdf(stringp);
end;
end



        
