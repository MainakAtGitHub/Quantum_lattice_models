function [deltaCal]= delta_full(eVector,fermi,nBands,nOrb,N,int,latt,deltaCal)
% old code        deltaCal = SCInteractionMatrix.*((eVector(1:nBands,:)*(((eVector((nBands + 1):end,:))').*repmat(fermi,1,nBands))));
% new code with for loops
% indices for deltaCal
parfor mu=1:nBands
    szlatt=size(latt);
    tmpvar=zeros(nOrb,1)+1i*zeros(nOrb,1);
    d=int32([0,0]);
    ymu=int32(mod(idivide(mu-1,int32(nOrb)),N)+1);
    xmu=int32(idivide(mu-1,int32(N*nOrb))+1);
    l1=mod(mu-1,nOrb)+1;
    for nu=1:nBands
        ynu=int32(mod(idivide(nu-1,int32(nOrb)),N)+1);
        xnu=int32(idivide(nu-1,int32(N*nOrb))+1);
        d(1)=xmu-xnu;
        d(2)=ymu-ynu;        % indices for bands
        if d(1)>= N/2
            d(1)=d(1)-N;
        end
        if d(2)>= N/2
            d(2)=d(2)-N;
        end        
        if d(1)<= -N/2
            d(1)=d(1)+N;
        end
        if d(2)<= -N/2
            d(2)=d(2)+N;
        end
        ind=int32(find(sum(latt==repmat(d,szlatt(1),1),2)==2));
        %if nu==85
         %   nu
        %end
        if ~isempty(ind)
            mu2=(ymu-1)*nOrb+(xmu-1)*nOrb*N;
            nu3=(ynu-1)*nOrb+(xnu-1)*nOrb*N;
            l4=mod(nu-1,nOrb)+1;            % indices for orbitals
            index0=(ind(1)-1)*nOrb^4+l1+(l4-1)*nOrb^3;
            %l2=l1;
            %l3=l4;
            for l2=1:nOrb
                % implement parfor here to get a speedup of a factor of
                % nOrb
                % remove the vectorization to get another speedup ?
                for l3 = 1:nOrb
                %for l3=1:nOrb
                    index=index0+(l2-1)*nOrb+(l3-1)*nOrb^2;
                    tmpvar(l3)=eVector(mu2+l2,:)*(eVector(nu3+l3+nBands,:)'.*fermi)*int(index);
                    %for n=1:nBands*2
                           %  SCInteractionMatrix(1,1).*((eVector(1,:)*(((eVector((nBands + 1),:))').*repmat(fermi,1,1))));
                        %deltaCal(mu,nu)=deltaCal(mu,nu)+eVector(mu2+l2,:)*(eVector(nu3+l3+nBands,:)'.*fermi)*int(index);%*eVector(mu2,n)*eVector(nu3+nBands,n)'*fermi(n);
                       % int(index)
                    %end;
                end
                deltaCal(mu,nu)=deltaCal(mu,nu)+sum(tmpvar);
            end
        end
    end
end







