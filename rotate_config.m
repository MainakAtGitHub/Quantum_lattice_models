function [rotated_config] = rotate_config(BdGfileName,uvec,th)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[filepath,name,ext]=fileparts(BdGfileName);
uvec=uvec/norm(uvec);
load(BdGfileName,'-mat');
if ~exist('nAnoUpDown','var')
    nAnoUpDown=zeros(size(nUp));
    nAnoDownUp=conj(nAnoUpDown);
end
save([BdGfileName,'[copied aside for rotation]'],'nAcc','nUpAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','energiesAcc','muAcc','mu','mudown','nUp','nDown','nUpdown','nDowndown','nAnoUpDown','nAnoDownUp');
magnetization=[(nAnoUpDown+nAnoDownUp),(1i*(nAnoUpDown-nAnoDownUp)),(nUp-nDown)]';
rot_mat=[cos(th)+uvec(1)^2*(1-cos(th)), uvec(1)*uvec(2)*(1-cos(th))-uvec(3)*sin(th), uvec(1)*uvec(3)*(1-cos(th))+uvec(2)*sin(th);...
    uvec(1)*uvec(2)*(1-cos(th))+uvec(3)*sin(th), cos(th)+uvec(2)^2*(1-cos(th)), uvec(2)*uvec(3)*(1-cos(th))-uvec(1)*sin(th);...
    uvec(1)*uvec(3)*(1-cos(th))-uvec(2)*sin(th), uvec(2)*uvec(3)*(1-cos(th))+uvec(1)*sin(th), cos(th)+uvec(3)^2*(1-cos(th))];
rotated_mag=(rot_mat*magnetization)';
nAnoUpDown=(rotated_mag(:,1)+rotated_mag(:,2)/(1i))/2;
nAnoDownUp=(rotated_mag(:,1)-rotated_mag(:,2)/(1i))/2;
nUp_new=((nUp+nDown)+rotated_mag(:,3))/2;
nDown_new=((nUp+nDown)-rotated_mag(:,3))/2;
nUp=nUp_new;
nDown=nDown_new;
save(BdGfileName,'nAcc','nUpAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','energiesAcc','muAcc','mu','mudown','nUp','nDown','nUpdown','nDowndown','nAnoUpDown','nAnoDownUp');
end

