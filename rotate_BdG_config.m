function f=rotate_BdG_config(inputfile)

read_input_file=inputfile;
read_input;
read_input_file

load(BdGfileName,'-mat');
rotated_nUp = zeros(size(nUp));
rotated_nDown = zeros(size(nDown));
for orbNo=1:nOrbitals
    sqr_nUp_for_rot = reshape(nUp(orbNo:nOrbitals:nOrbitals*N^2),[N,N]);
    sqr_nUp_for_rot = sqr_nUp_for_rot';
    rotated_nUp(orbNo:nOrbitals:nOrbitals*N^2) = sqr_nUp_for_rot(:);
    
    sqr_nDown_for_rot = reshape(nDown(orbNo:nOrbitals:nOrbitals*N^2),[N,N]);
    sqr_nDown_for_rot = sqr_nDown_for_rot';
    rotated_nDown(orbNo:nOrbitals:nOrbitals*N^2) = sqr_nDown_for_rot(:);
end
save([Gamma_file,'_rotated_for_n'],'delta','mu','nUp','nDown');