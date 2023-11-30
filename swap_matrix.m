function sw=swap_matrix(input)
if isa(input,'char')
    load(input);
    matrix=Impparameters;
else
    matrix=input;
end;
number=size(matrix,3);
sw=0*matrix;
nOrb=size(matrix,2);
nOrb2=nOrb/2;
if round(nOrb2*2)==nOrb
    for n=1:number
        sw(:,:,n)=matrix([(nOrb2+1):nOrb,1:nOrb2],[(nOrb2+1):nOrb,1:nOrb2],n);
    end;
else
    disp('odd number of orbitals, cannot swap sublattices')
end;
if isa(input,'char')
    newfile=[input,'_swap'];
    copyfile(input,newfile);
    m = matfile(newfile);
    m.Properties.Writable=true;
    m.Impparameters=sw;
end;
    