function t1=triD(e,k,E)
% e input energy for which value of density of states is desired
% k matrix of three k-vectors of corners of triangle
% E corresponding energies at the corners of triangle
%------------------------------------------------------
%sort energies in ascending order to simplify the decission process
A=sortrows([e,k]);
e=A(:,1);
k=A(:,2:size(k,2)+1);
%define difference vectors
k12=k(2,:)-k(1,:);
k13=k(3,:)-k(1,:);
k23=k(3,:)-k(2,:);
%calculate the gradient on the triangle (is constant over triangle)
 if e(2)==e(1)&&e(3)==e(2)&&e(3)==e(1)
     %avoid numerical problems with flat dispersion
     gradE=eps
 else
     gradE=[k12(2)*(e(3)-e(1))-(e(2)-e(1))*k13(2),(e(2)-e(1))*k13(1)-k12(1)*(e(3)-e(1))]/(k12(1)*k13(2)-k12(2)*k13(1));
 end;
% elseif e(2)==e(3)
%     gradE=k12*(e(2)-e(1))/norm(k12)^2+k13*(e(3)-e(1))/norm(k13)^2;
% elseif e(2)==e(1)
%     gradE=k13*(e(3)-e(1))/norm(k13)^2+k23*(e(3)-e(2))/norm(k23)^2;
% else
%     gradE=k12*(e(2)-e(1))/norm(k12)^2+k13*(e(3)-e(1))/norm(k13)^2;
% end;
%calculate absolute value of gradient
de=norm(gradE);
%initialize the result vector
%t1=zeros(1,length(E));
%calculate dos for every engergy point with prepared input arguments
%for i=1:length(E)
%    t1(i)=tri(E(i),k12,k13,k23,e)/de;
%end;
% vectorized version
t1=triV(E,k12,k13,k23,e)/de;
%if sum(t1(:)>0)
%if ~isequal(t1,t1E)
%    t1(1)
%end;
%end;

%calculation of dos (without normalization) for one energy point
function t=tri(E,k12,k13,k23,e)
if E <= e(1)
    %energy below lowest energy: no states avalable
    t=0;
elseif E <= e(2)
    % density of states proportional to length of cutting line in triangle
    % according the forumula similar to Bloechel94
    a=(E-e(1))/(e(2)-e(1))*k12;
    b=(E-e(1))/(e(3)-e(1))*k13;
    t=norm(a-b);
elseif E <= e(3)
    c=(E-e(2))/(e(3)-e(2))*k23;
    b=(E-e(1))/(e(3)-e(1))*k13;
    t=norm(k12+c-b);
else
    %energy above highest energy of triangle: no states
    t=0;
end;
end;

function t=triV(E,k12,k13,k23,e)
E_smaller_e2=(E > e(1)).*(E<=e(2));
E_smaller_e3=(E>e(2)).*(E<=e(3));
lE=length(E);
t=zeros(1,lE);
%if E <= e(1)
    %energy below lowest energy: no states avalable
%    t=0;
%elseif E <= e(2)
    % density of states proportional to length of cutting line in triangle
    % according the forumula similar to Bloechel94
    a=k12'*(E_smaller_e2.*(E-e(1))/(e(2)-e(1)));
    b=k13'*(E_smaller_e2.*(E-e(1))/(e(3)-e(1)));
    t=t+sqrt(sum((a-b).^2,1));
%elseif E <= e(3)
    c=k23'*(E_smaller_e3.*(E-e(2))/(e(3)-e(2)));
    b=k13'*(E_smaller_e3.*(E-e(1))/(e(3)-e(1)));
    t=t+sqrt(sum((repmat(k12',1,lE).*repmat(E_smaller_e3,2,1)+c-b).^2,1));
%else
    %energy above highest energy of triangle: no states
%    t=0;
end
end