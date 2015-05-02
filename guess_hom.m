function r=guess_hom(inputfile,newN)
        % new text-based input format
        read_input_file=inputfile;
        read_input;
        read_input_file
        if nargin <2
            newN=input('new size: ');
        end;
        global fullgamma
% spin polarized calculation
if ~exist('spinpolarized','var')
    spinpolarized=false;
end;     
load(BdGfileName,'-mat');
load(Gamma_file,'-mat')
if exist('Gamma','var')
    fullgamma=false;
else
    fullgamma=true;
    latticeVectorsSC=Gammafull.latt;
end;
delta_hom=homogenize_delta(delta,latticeVectorsSC,nOrb,N,newN);
nup1=nUp(1:nOrbitals);
ndown1=nDown(1:nOrbitals);
cutoff=1e-17;
cutoff=eps;
if exist('nup1','var')
    nUp=repmat(nup1,newN^2,1);
    nDown=repmat(ndown1,newN^2,1);
    % clean values
   delta=(abs(real(delta_hom))>cutoff).*real(delta_hom)+(abs(imag(delta_hom))>cutoff).*imag(delta_hom);
end;
 deltaMaxAcc = [];
 deltaMinAcc = [];
 deltaDiffAcc = [];
 muAcc = [];
 nAcc=[];
    if ~spinpolarized
        save([BdGfileName,num2str(newN)],'nAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','muAcc','mu', 'nUp','nDown');
    else
        save([BdGfileName,num2str(newN)],'nAcc','delta','deltaMaxAcc','deltaMinAcc','deltaDiffAcc','muAcc','mu','mudown','nUp','nDown','nUpdown','nDowndown');
    end