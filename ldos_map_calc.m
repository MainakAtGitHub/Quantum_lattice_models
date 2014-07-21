function localLdos=ldos_map_calc(N,nOrbitals,xGridRange_limits,yGridRange_limits,zGridPoint,RDiscrete,shift,sizeWannier,wannierValuesreshape,latticeGreens)
minarg=sizeWannier(1)*sizeWannier(2)*sizeWannier(3)
maxarg=0
numlm=N^2;%(2*(ceil(N/2)-1)+1)^2;
%wAcc=zeros(numlm,nOrbitals);
lhalf=ceil(N/2)-1;
mhalf=ceil(N/2)-1;
xGridRange=xGridRange_limits(1):xGridRange_limits(2);
yGridRange=yGridRange_limits(1):yGridRange_limits(2);
% calculation part of ldos_map
    localLdos = zeros(length(xGridRange),length(yGridRange));
    countLoopX = 0;
    for xGridPoint = xGridRange_limits(1):xGridRange_limits(2);
        disp(['Calculating ',my_int2str(xGridPoint), 'of (',my_int2str(xGridRange(1)),'..',my_int2str(xGridRange(numel(xGridRange))),')']);
        countLoopX = countLoopX + 1;
        countLoopY = 0;
        for yGridPoint =yGridRange_limits(1):yGridRange_limits(2);
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
            [l1,m1]=meshgrid(lrange,mrange);
            l1=l1(:);
            m1=m1(:);
            szl=numel(l1);
            %position=N*(mod(l1+lhalf,N))+mod(m1+mhalf,N)+1;
            position=(mod(l1+lhalf,N))+N*mod(m1+mhalf,N)+1;
            shiftr=repmat(shift,szl,1);
            rr=repmat(r,szl,1);
            RDiscreter=repmat(RDiscrete,szl,1);
            shiftedArgument = rr -  RDiscreter.*[l1 m1 zeros(szl,1)] + shiftr;
            indexshiftedArgument=shiftedArgument(:,1)+(shiftedArgument(:,2)-1)*sizeWannier(1)+(shiftedArgument(:,3)-1)*sizeWannier(1)*sizeWannier(2);
            minarg= min(minarg,min(indexshiftedArgument(:)));
            maxarg= max(maxarg,max(indexshiftedArgument(:)));
            % for n=1:szl
            % remove loop totally
            %for l = lrange
            %   for m = mrange
            %for l = -(ceil(N/2)-1):(ceil(N/2)-1)%lrange
            %   for m = -(ceil(N/2)-1):(ceil(N/2)-1)
            %position=position+1;
            %position=N*(mod(l1(n)+lhalf,N))+mod(m1(n)+mhalf,N)+1;
            %if position>numlm
            %    disp(['outside']);
            %end;
            %R = [l m 0];
            %latticeVector = RDiscrete.*R;
            %wannierArgument = r - latticeVector;
            %shiftedArgument = wannierArgument + shift; % translate wannier origin
            % shift to central area if needed (for x and y direction)
            % check whether this argument is in range or not
            %if (((shiftedArgument) <= sizeWannier) & ((shiftedArgument) >= [1 1 1]))
            % yes in range, now find the value
            %wannierValue = wannierI(shiftedArgument);
            %for orbital = 1:nOrbitals
            %                 w = squeeze(wannierValues(shiftedArgument(n,1),shiftedArgument(n,2),shiftedArgument(n,3),:));
            %  w1=squeeze(wannierValuesreshape(:,indexshiftedArgument));
            wAcc(:,position) = wannierValuesreshape(:,indexshiftedArgument);
            %end;
            %else
            %disp(['outside:',num2str(shiftedArgument)]);
            % not in range, set it to zero
            % wAcc = [wAcc; zeros(nOrbitals,1)];
            % end
            %end
            %end
            %  wAcc=wAcc';
            wAcc1=wAcc(:);
            localLdos(countLoopX,countLoopY) = (-1/pi)*imag(wAcc1'*(latticeGreens*wAcc1));
        end
    end
    minarg
    maxarg
