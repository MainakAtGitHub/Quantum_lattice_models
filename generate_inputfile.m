%generate inputfiles for various kT, and U,U'=U/2,J=U/4,J'=0 (for now, 8March2020)
function f=generate_inputfile(oldFolder, Gamma_file,... 
    n0_orb,random_n, neel_n, stripe_n,...
    pos_file, ...
    N, Vimp,kT, U, U_Upr_J_related,Upr, J,...
    nOrbitals,TB_file,nTol,justRun,maxLoop,...
    ref_grid_hopping_file)
% if any(~exist('random_n','var'),random_n==[],nargin<4)
%     random_n=false;
% end
% if any(~exist('neel_n','var'),neel_n==[],nargin<5)
%     neel_n=false;
% end
% if any(~exist('stripe_n','var'),stripe_n==[],nargin<6)
%     stripe_n=false;
% end
% if any(~exist('Vimp','var'),Vimp==[],nargin<9)
%     Vimp=0;
% end
% if any(~exist('U_Upr_J_related','var'),U_Upr_J_related==[],nargin<12)
%     U_Upr_J_related=false;
% end
% if any(~exist('Upr','var'),Upr==[],nargin<13)
%     Upr=0;
% end
% if any(~exist('J','var'),J==[],nargin<14)
%     J=0;
% end
    
BdGfolder='/home/UFAD/mainak.pal/Desktop/summer_2019_office/andreas/BdG';


% kT = 0.06:0.02:0.1;
% U = 2:0.1:2.5;
for itr_N=N
    for itr_Vimp=Vimp
        for itr_kT = kT
            for itr_U = U
                if U_Upr_J_related
                    Upr=itr_U/2;
                end
                if U_Upr_J_related
                    J=itr_U/4;
                end
                cd(oldFolder);
                inputfile = [...
                    'input_1band_N',num2str(itr_N),...
                    '_random_initial_strength',num2str(random_n),...
                    '_neel_initial_strength',num2str(neel_n),...
                    '_stripe_initial_strength',num2str(stripe_n),...
                    '_Vimp',num2str(itr_Vimp),...
                    '_U',num2str(itr_U),...
                    '_Upr',num2str(Upr),...
                    '_J',num2str(J),...
                    '_kT',num2str(itr_kT),...
                    '_n0',num2str(n0_orb*nOrbitals),'.txt'];
                infile_id = fopen(inputfile,'w');
                
                BdGfilename=[...
                    'BdG_1band_N',num2str(itr_N),...
                    '_random_initial_strength',num2str(random_n),...
                    '_neel_initial_strength',num2str(neel_n),...
                    '_stripe_initial_strength',num2str(stripe_n),...
                    '_Vimp',num2str(itr_Vimp),...
                    '_U',num2str(itr_U),...
                    '_Upr',num2str(Upr),...
                    '_J',num2str(J),...
                    '_kT',num2str(itr_kT),...
                    '_n0',num2str(n0_orb*nOrbitals),'.mat'];
                fprintf(infile_id,['BdGfileName=',num2str(BdGfilename),'\n']);
                [Gamma_filepath,Gamma_filename,Gamma_fileext]=fileparts(Gamma_file);                     
                fprintf(infile_id,['Gamma_file=',num2str(Gamma_filename),'.mat','\n']);
                fprintf(infile_id,'M=5\n');
                fprintf(infile_id,['N=',num2str(itr_N),'\n']);
                [TB_filepath,TB_filename,TB_fileext]=fileparts(TB_file);
                fprintf(infile_id,['TB_file=',num2str(TB_filename),'.mat','\n']);
                fprintf(infile_id,['Vimp=',num2str(itr_Vimp),'\n']);   % potential loop
                fprintf(infile_id,'alpha=0.25\n');
                fprintf(infile_id,'beta1=0.7\n');
                fprintf(infile_id,'beta2=0.7\n');
                fprintf(infile_id,'casename=Mainak\n');
                fprintf(infile_id,['casestring=','LDOS_',...
                    'input_1band_N',num2str(N),...
                    '_random_initial_strength',num2str(random_n),...
                    '_neel_initial_strength',num2str(neel_n),...
                    '_stripe_initial_strength',num2str(stripe_n),...
                    '_Vimp',num2str(itr_Vimp),...
                    '_U',num2str(itr_U),...
                    '_Upr',num2str(Upr),...
                    '_J',num2str(J),...
                    '_kT',num2str(itr_kT),...
                    '_n0',num2str(n0_orb*nOrbitals),'\n']);
                fprintf(infile_id,'deltaTol=1e-06\n'); % may need frequent change
                fprintf(infile_id,'firstEnergy=-5\n');
                fprintf(infile_id,'ita=0.005\n');
                fprintf(infile_id,['kT=',num2str(itr_kT),'\n']);
                fprintf(infile_id,'lastEnergy=5\n');
                fprintf(infile_id,['maxLoop=',num2str(maxLoop),'\n']);
                fprintf(infile_id,['n0=',num2str(n0_orb*nOrbitals),'\n']);
                fprintf(infile_id,'nEnergyPoints=2500\n');
                fprintf(infile_id,['nOrbitals=',num2str(nOrbitals),'\n']);
                fprintf(infile_id,'tetra=false\n');
                fprintf(infile_id,'calcGreens=false\n');
                fprintf(infile_id,'Greensenergy=-0.0036\n');
                fprintf(infile_id,['nTol=',num2str(nTol),'\n']);
                fprintf(infile_id,'singular_quad=true\n');
                fprintf(infile_id,'correlated=true\n');
                fprintf(infile_id,'spinpolarized=true\n');
                fprintf(infile_id,['U=',num2str(itr_U),'\n']);
                fprintf(infile_id,['U_pr=',num2str(Upr),'\n']);
                fprintf(infile_id,['J=',num2str(J),'\n']);
                fprintf(infile_id,'saveEigVec=false\n');
                fprintf(infile_id,'correlated=true\n');
                
                if true;%and(exist('pos_file','var'),pos_file~=[])
                    [pos_filepath,pos_filename,pos_fileext]=fileparts(pos_file);
                    
                    fprintf(infile_id,['pos_file=',...
                        num2str(pos_filename),num2str(pos_fileext),'\n']);
                end
                fprintf(infile_id,'normal_metal=true\n');
                
                [ref_grid_hopping_filepath,ref_grid_hopping_filename,ref_grid_hopping_fileext]=fileparts(ref_grid_hopping_file);
                fprintf(infile_id,['ref_grid_hopping_file=',...
                    num2str(ref_grid_hopping_filename),num2str(ref_grid_hopping_fileext),'\n']);
                
                fclose(infile_id);
                if ~justRun
                    submit_inputfile = [...
                        'submit_input_1band_N',num2str(itr_N),...
                        '_random_initial_strength',num2str(random_n),...
                        '_neel_initial_strength',num2str(neel_n),...
                        '_stripe_initial_strength',num2str(stripe_n),...
                        '_Vimp',num2str(itr_Vimp),...
                        '_U',num2str(itr_U),...
                        '_Upr',num2str(Upr),...
                        '_J',num2str(J),...
                        '_kT',num2str(itr_kT),...
                        '_n0',num2str(n0_orb*nOrbitals),'.sh'];
                    submit_infile_id = fopen(submit_inputfile,'w');
                    fprintf(submit_infile_id,'#!/bin/bash\n');
                    fprintf(submit_infile_id,'#SBATCH --job-name=parallel_job      # Job name\n');
                    fprintf(submit_infile_id,'##SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)\n');
                    fprintf(submit_infile_id,'##SBATCH --mail-user=email@ufl.edu    # Where to send mail\n');
                    fprintf(submit_infile_id,'#SBATCH --ntasks=1                   # Run a single task\n');
                    fprintf(submit_infile_id,'#SBATCH --cpus-per-task=4            # Number of CPU cores per task\n');
                    fprintf(submit_infile_id,'#SBATCH --mem=6gb                    # Job memory request\n');
                    fprintf(submit_infile_id,'#SBATCH --time=10:00:00              # Time limit hrs:min:sec\n');
                    fprintf(submit_infile_id,'#SBATCH --output=slurm_%N_%j.out     # Standard output and error log\n');
                    fprintf(submit_infile_id,'##SBATCH --partition=hpg2-dev\n');
                    fprintf(submit_infile_id,['./run_BdG_impurity_v3.sh /apps/matlab/r2019b ',inputfile]);
                    fclose(submit_infile_id);
                end
                
                % now copy initialguess to bdgfile, *****IF adjusting nTol MUST
                % MUST MUST
                % comment LINE BELOW
%               cd(oldFolder);
                if ~justRun
                    cd(BdGfolder);
                    generate_initialguess(itr_N,Gamma_file,n0_orb,0,...
                        random_n,neel_n,stripe_n); %mu = 0 initial guess
                    cd(oldFolder);
                    ini_guess_filename=[Gamma_file,'guess',num2str(itr_N),...
                        '_rand_strength_',num2str(random_n),...
                        '_neel_strengh_',num2str(neel_n),...
                        '_stripe_strengh_',num2str(stripe_n)];
                    copyfile(ini_guess_filename,BdGfilename);
                    copyfile(ini_guess_filename,['Ini_config_',BdGfilename]);
                end
%               copyfile('Gamma_0_NN.matguess12_rand_strength_0_neel_strengh_0.01_stripe_strengh_0',...
% ['input_5band_N12_neel_initial_Vimp0_U',num2str(itr_U),'_Upr',num2str(itr_U/2),'_J',num2str(itr_U/4),...
% '_kT',num2str(itr_kT),'_n0_6.mat'])
%                 if false;%or(~exist('pos_file','var'),pos_file==[])
%                     dislocation_length=0;
%                 else
                    load(pos_file,'-mat');
                    dislocation_length=itr_N^2-size(r,1);
                    if dislocation_length>0
                        load(BdGfilename,'-mat');
                        delta=zeros(size(r,1));
                        nUp([size(r,1)+1:end])=[];
                        nDown([size(r,1)+1:end])=[];
                        save(BdGfilename,'delta','nUp','nDown','mu');
                    end
                    if dislocation_length>0
                        load(['Ini_config_',BdGfilename],'-mat');
                        delta=zeros(size(r,1));
                        nUp([size(r,1)+1:end])=[];
                        nDown([size(r,1)+1:end])=[];
                        save(['Ini_config_',BdGfilename],'delta','nUp','nDown','mu');
                    end
                    
%                 end
                  cd(BdGfolder);
                  inputfileWithPath=[oldFolder,'/',inputfile];
                  if justRun
                      BdG_impurity_v3(inputfileWithPath,0,false);
                  end
                  close all;
                  plot_mag_den_1band(inputfileWithPath);
            end
        end
    end
end
cd(BdGfolder);
end
