function f1=triF(e,a,k,E)
% E input energy for which value of density of states is desired
% a function to be integrated (given at the edge points of the triangle)
% k matrix of three k-vectors of corners of triangle
% e corresponding energies at the corners of triangle
%------------------------------------------------------
%sort energies in ascending order to simplify the decission process
[e,index]=sortrows(e);
%e1=B;
a=a(index,:);
k=k(index,:);
%A=sortrows([e,a,k]);
%e=A(:,1);
%a=A(:,2);
%k=A(:,3:size(k,2)+2);
%define difference vectors
k12=k(2,:)-k(1,:);
k13=k(3,:)-k(1,:);
k23=k(3,:)-k(2,:);
% calculate the coeficients to approximate the function a over the triangle
a0=a(1,:);%(k(1,1)*(k(2,2)*a(3)-a(2)*k(3,2))+k(1,2)*(a(2)*k(3,1)-k(2,1)*a(3))+a(1)*(k(2,1)*k(3,2)-k(2,2)*k(3,1)))/(k12(1)*k13(2)-k12(2)*k13(1));
a1=-(k12(2)*(a(3,:)-a(1,:))-(a(2,:)-a(1,:))*k13(2))/(k12(1)*k13(2)-k12(2)*k13(1));
% debug code
if sum(isnan(a1)>0)
    disp('a0 isnan');
end;
a2=-((a(2,:)-a(1,:))*k13(1)-k12(1)*(a(3,:)-a(1,:)))/(k12(1)*k13(2)-k12(2)*k13(1));
% debug code
if sum(isnan(a2)>0)
    disp('a2 isnan');
end;
% calculate the density of states
i0=triD(e,k,E);
% calculate the momentum integrals distinguishing the tree cases
i=triI(e,k,E);
f1=a0'*i0+a1'*i(1,:)+a2'*i(2,:);


function t1=triI(e,k,E)
%define difference vectors
%k12=k(2,:)-k(1,:);
%k13=k(3,:)-k(1,:);
%k23=k(3,:)-k(2,:);
%calculate the gradient on the triangle (is constant over triangle)
 if e(2)==e(1)&&e(3)==e(2)&&e(3)==e(1)
     %avoid numerical problems with flat dispersion
     gradE=eps
 else
     gradE=[k12(2)*(e(3)-e(1))-(e(2)-e(1))*k13(2),(e(2)-e(1))*k13(1)-k12(1)*(e(3)-e(1))]/(k12(1)*k13(2)-k12(2)*k13(1));
 end;
%calculate absolute value of gradient
 de=norm(gradE);
%initialize the result vector
%t1=zeros(2,length(E));
%calculate dos for every engergy point with prepared input arguments
%for i=1:length(E)
%    t1(:,i)=tri(E(i),k12,k13,k23,e)/de;
%end;
% vectorized version
t1=triV(E,k12,k13,k23,e)/de;
%if sum(sum(t1(:)))>0
%if ~isequal(t1,t1E)
%    t1(1)
%end;
%end;
end;


function t=tri(E,k12,k13,k23,e)
if E <= e(1)
    %energy below lowest energy: no states avalable
    t=[0;0];
elseif E <= e(2)
    % density of states proportional to length of cutting line in triangle
    % according the forumula similar to Bloechel94
    av=(E-e(1))/(e(2)-e(1))*k12;
    b=(E-e(1))/(e(3)-e(1))*k13;
    t=norm(av-b)*1/2*(av+b)';
elseif E < e(3)
    c=(E-e(2))/(e(3)-e(2))*k23;
    b=(E-e(1))/(e(3)-e(1))*k13;
    t=1/2*norm(k12+c-b)*(k12+c+b)';
else
    %energy above highest energy of triangle: no states
    t=[0;0];
end;
end

function t=triV(E,k12,k13,k23,e)
lE=length(E);
t=zeros(2,lE);
%a=zeros(1,lE);
%b=zeros(1,lE);
%k12=repmat(k12,lE);
%k13=repmat(k13,lE);
%k23=repmat(k23,lE);
%if E <= e(1)
    %energy below lowest energy: no states avalable
%    t=[0;0];
E_smaller_e2=(E > e(1)).*(E<=e(2));
E_smaller_e3=(E>e(2)).*(E<=e(3));
%E_bigger_e1=(E > e(1));
%E_smaller_e2=E_bigger_e1.*(E<=e(2));
%elseif E <= e(2)
    % density of states proportional to length of cutting line in triangle
    % according the forumula similar to Bloechel94
    av=k12'*(E_smaller_e2.*(E-e(1))/(e(2)-e(1)));
%    a=(E-e(1))/(e(2)-e(1))*k12;
    b=k13'*(E_smaller_e2.*(E-e(1))/(e(3)-e(1)));
    t=t+repmat(sqrt(sum((av-b).^2,1)),2,1).*1/2.*(av+b);
    av=k23'*(E_smaller_e3.*(E-e(2))/(e(3)-e(2)));
    b=k13'*(E_smaller_e3.*(E-e(1))/(e(3)-e(1)));
    t=t+1/2*repmat(sqrt(sum((repmat(k12',1,lE).*repmat(E_smaller_e3,2,1)+av-b).^2,1)),2,1).*(repmat(k12',1,lE)+av+b);
    % remove nan numbers that come from 0*1/(0) if e(1)==e(2) etc.
    t(isnan(t))=0;
%else
    %energy above highest energy of triangle: no states
%    t=[0;0];
end
end