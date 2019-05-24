function h=colorbar_rwb(figure1,maxabsekkn,ticks,labels,abs)
global colorred
    global fsz;
 if isempty(fsz)
 fsz=10;
 end;
color1=[1 0 0]; % red
%color2=[0 0 0]; % black
color2=[1 1 1]; % white
color3=[0 0 1]; % blue
if ~exist('colorred','var')
    colorrred=false;
end;
if ~colorred
    input=[-10 -1 0 1 10]+.5;
    colormatrix=[color1; color1;color2;color3; color3];
else
    input=[0 1 10];
    colormatrix=[color2;color1; color1];
end;
if ~isempty(figure1)
    % draw a colorbar
%     if ~colorred
%     x=-maxabsekkn:maxabsekkn/255:maxabsekkn;
%     else
            x=0:maxabsekkn/255:maxabsekkn;
%    end
    r=interp1(input*maxabsekkn,colormatrix(:,1),x);
g=interp1(input*maxabsekkn,colormatrix(:,2),x);
b=interp1(input*maxabsekkn,colormatrix(:,3),x);
    set(figure1,'Colormap', [r',g',b']);
        allAxesInFigure = findall(figure1,'type','axes');
if nargin >2
    h = colorbar;
    if exist('abs','var')
        ticks_res=ticks;
    else
        ticks_res=round(ticks/maxabsekkn*128+128);
    end
    if nargin < 4
        labels = num2str(repmat(ticks, 1, 1)', 2);
    end
    % eliminate the same ticks_res
    ticksres1=ticks_res(1);
    labels1=labels(1,:);
    for n=2:length(ticks_res)
        if ticks_res(n)> ticksres1(length(ticksres1))
            ticksres1=[ticksres1,ticks_res(n)];
            labels1=[labels1;labels(n,:)];
        end
    end;
    set(h, 'YTick', ticksres1);
set(h, 'YTickLabel', labels1);
if ~colorred
    set(allAxesInFigure,'FontSize',fsz); 
else
        set(allAxesInFigure,'FontSize',fsz); 
end;
else

if ~colorred
    set(allAxesInFigure,'CLim',[-maxabsekkn maxabsekkn],'FontSize',fsz); 
else
        set(allAxesInFigure,'CLim',[0 maxabsekkn],'FontSize',fsz); 
end;
end

end;

   % cb=colorbar('peer',allAxesInFigure,'FontSize',fsz);