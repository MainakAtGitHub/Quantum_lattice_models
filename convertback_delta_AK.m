function c=convertback_delta_AK(filename,N,nOrb,fill,mu)
if nargin <2
    N=9;
end
if nargin <3
    nOrb=10;
end;
if nargin < 4
    fill=0.6;
end;
if nargin < 5
    mu=0.03;
end;
r=read_complex_matrix(filename,nOrb,nOrb);
% we assume to have the 11x11 spatial grid
delta=zeros(N^2*nOrb);
for ix=1:N
    for iy=1:N
        for jx=1:N
            for jy=1:N
                indexx=ix-jx+6;
                indexy=iy-jy+6;
                if (indexx>0) && (indexy>0) && (indexx<12) && (indexy<12)                    
                    index=indexx+11*(indexy-1);
                    [iRange, jRange]=find_lattice_translation_index(N,nOrb,[ix iy],[jx jy]);
                    delta(iRange,jRange)=reshape(r(:,index),nOrb,nOrb);
                end;
                
            end
        end
    end
end
% initial guess for the densities
nUp=ones(N^2*nOrb,1)*fill;
nDown=ones(N^2*nOrb,1)*fill;
save([filename,'convb'],'delta');