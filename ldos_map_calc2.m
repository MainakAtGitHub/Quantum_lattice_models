function localLdos=ldos_map_calc2(N,nOrbitals,xGridRange_limits,yGridRange_limits,zGridPoint,RDiscrete,shift,sizeWannier,wannierValuesreshape,latticeGreens)
numlm=N^2;%(2*(ceil(N/2)-1)+1)^2;
%wAcc=zeros(numlm,nOrbitals);
lhalf=fix(N/2);
mhalf=fix(N/2);
xGridRange=xGridRange_limits(1):xGridRange_limits(2);
yGridRange=yGridRange_limits(1):yGridRange_limits(2);
    localLdos = zeros(length(xGridRange),length(yGridRange));
    countLoopX = 0;
    for xGridPoint = xGridRange
        disp(['Calculating ',my_int2str(xGridPoint), 'of (',my_int2str(xGridRange(1)),'..',my_int2str(xGridRange(numel(xGridRange))),')']);
        countLoopX = countLoopX + 1;
        countLoopY = 0;
        for yGridPoint = yGridRange
            countLoopY = countLoopY + 1;
            r = [xGridPoint, yGridPoint, zGridPoint];
           % wAcc = [];
           %position=0;
           wAcc=zeros(nOrbitals,numlm);
           % change way of calculation to reduce loop
           % calculation of minimal / maximal l and m where we have Wannier
           % boxes
           centerm=yGridPoint/RDiscrete(1);
           centerl=xGridPoint/RDiscrete(1);
           minl=round(centerl-(sizeWannier(1)-shift(1)+0.5)/RDiscrete(1)+0.5);
           maxl=round(centerl+(sizeWannier(1)-shift(1)+0.5)/RDiscrete(1)-0.5);
           minm=round(centerm-(sizeWannier(2)-shift(2)+0.5)/RDiscrete(2)+0.5);
           maxm=round(centerm+(sizeWannier(2)-shift(2)+0.5)/RDiscrete(2)-0.5);
           lrange=(minl:maxl);
           mrange=(minm:maxm);
%            for l = lrange
%                for m = mrange
%                    %for l = -(ceil(N/2)-1):(ceil(N/2)-1)%lrange
%                    %   for m = -(ceil(N/2)-1):(ceil(N/2)-1)
%                    %position=position+1;
%                    position=N*(mod(l+lhalf,N))+mod(m+mhalf,N)+1;
%                    %if position>numlm
%                    %    disp(['outside']);
%                    %end;
%                    %R = [l m 0];
%                    %latticeVector = RDiscrete.*R;
%                    %wannierArgument = r - latticeVector;
%                    %shiftedArgument = wannierArgument + shift; % translate wannier origin
%                    % shift to central area if needed (for x and y direction)
%                    shiftedArgument = r -  RDiscrete.*[l m 0] + shift;
%                    % check whether this argument is in range or not
%                    %if (((shiftedArgument) <= sizeWannier) & ((shiftedArgument) >= [1 1 1]))
%                    % yes in range, now find the value
%                    %wannierValue = wannierI(shiftedArgument);
%                    %for orbital = 1:nOrbitals
%                    indexshiftedArgument=shiftedArgument(:,1)+(shiftedArgument(:,2)-1)*sizeWannier(1)+(shiftedArgument(:,3)-1)*sizeWannier(1)*sizeWannier(2);
%                    %                            w = squeeze(wannierValuesreshape(shiftedArgument(1),shiftedArgument(2),shiftedArgument(3),:));
%                    wAcc(:,position) = wannierValuesreshape(:,indexshiftedArgument);
%                    %                           wAcc(position,:) =  w;
%                    %end;
%                    %else
%                    %disp(['outside:',num2str(shiftedArgument)]);
%                    % not in range, set it to zero
%                    % wAcc = [wAcc; zeros(nOrbitals,1)];
%                    % end
%                end
%            end
%            wAcc=wAcc(:);

% remove vectorization in order to take into account the
            % nonlocal contributions at the edges correctly
            % orbitals not yet implemented
            for l = lrange
                for m = mrange
                    xi=N*(mod(l+lhalf,N))+mod(m+mhalf,N)+1;
                    shiftedArgument = r -  RDiscrete.*[l m 0] + shift;
                    indexshiftedArgument=shiftedArgument(:,1)+(shiftedArgument(:,2)-1)*sizeWannier(1)+(shiftedArgument(:,3)-1)*sizeWannier(1)*sizeWannier(2);
                    w= wannierValuesreshape(:,indexshiftedArgument);
                    for l1 = lrange
                        for m1 = mrange
                            yi=N*(mod(l1+lhalf,N))+mod(m1+mhalf,N)+1;
                            shiftedArgument1 = r -  RDiscrete.*[l1 m1 0] + shift;
                            indexshiftedArgument1=shiftedArgument1(:,1)+(shiftedArgument1(:,2)-1)*sizeWannier(1)+(shiftedArgument1(:,3)-1)*sizeWannier(1)*sizeWannier(2);
                            w1= wannierValuesreshape(:,indexshiftedArgument1);
                            localLdos(countLoopX,countLoopY) =localLdos(countLoopX,countLoopY) -1/pi*imag(w*latticeGreens(xi,yi)*w1);
                        end
                    end
                end
            end
                                      
        %    for xi=1:numlm
         %       for yi=1:numlm
         %           localLdos(countLoopX,countLoopY) =localLdos(countLoopX,countLoopY) -1/pi*imag(wAcc(yi)*latticeGreens(xi,yi)*wAcc(xi));
         %       end
         %   end;
        end
    end