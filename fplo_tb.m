%% -- dx2mat.m --
% 
% this function reads in a DX file and translates it into a 3D matrix in
% matlab. DX files are file formats for storing 3D data in plaintext, and
% they are readable by PyMOL and VMD. One inputs the DX file path, and the
% program outputs a structure with all the information about the DX
% 
% As of July 2012, file output format is as defined on
% http://www.poissonboltzmann.org/file-formats/mesh-and-data-formats/opendx-scalar-data
% 
% The output of this can be read by mat2dx.m
% 
% -- required inputs (1) --
%
% input value           meaning
%
% file name             a string that leads to the DX file path
%
% -- outputs (5) --
%
% output value         meaning                           example value
%
% DXdata.densityMatrix  rectangular-prism 3D matrix holding density values.
%                       All XYZ values correspond to spatial location of
%                       data.
% 
% DXdata.outfile      output file name                 "mat2dx.dx"
% DXdata.minX         origin (min) X coord (angstroms)  0.00
% DXdata.minY         origin (min) Y coord (angstroms)  0.00
% DXdata.minZ         origin (min) Z coord (angstroms)  0.00
% DXdata.voxelLength  length to side of voxel (angst.)  1.00 or [1.10 2.30 ]
% 
% 
% -- example usage: normalize DX output --
% 
% DXdata = dx2mat('DXfile.dx');                      % interpret DX data
% mat3D = DXdata.densityMatrix;                      % make name convenient
% maxValue = max(max(max(mat3D)));                   % find max value
% minValue = min(min(min(mat3D)));                   % find min value
% valueRange = maxValue - minValue;                  % find range of data
% DXdata.densityMatrix = (mat3D-minValue)/valueRange; % normalize and save
% DXdata.outfile = 'DXfile_normaized.dx';            % rename output
% mat2dx(DXdata);

function [tb,latt_vec] = fplo_tb(outfile,number,complx)
tic;
if nargin < 2
    % give back last tb model
    number=-1;
end
if nargin < 3
    % give back last tb model
    complx=0;
end;
%% read every line in file

% initialize file (read all into memory)
fprintf('reading in file %s\n', outfile);
FileID = fopen(outfile);
rawText = fread(FileID,inf,'*char');

% close file
fclose(FileID);

%% interpret introductory data

fprintf('importing data about 3D matrix\n');

% parse lines by end-of-lines
splitLines = textscan(rawText, '%s', 'delimiter', '\n');
splitLines=splitLines{1};
%splitLines = strread(rawText, '%s', 'delimiter', '\n');
clear rawText
numLines = length(splitLines);
% find positions of tb models in output file

pos=false;
% use introductory data to deduce file length and grid positions
voxelLength = [0 0 0]';
lineCounter = 1;
% get rid of the loop by use of strncmp(string,cellstr,n)
% startIndex = regexp(str,expression)

% read in lattice constants
str0='lattice constants [aB] :';
lat_pos=strncmp(str0,splitLines,length(str0));
str1='lattice vectors';
lat_pos1=strncmp(str1,splitLines,length(str1));
lat_pos_index1=find(lat_pos1==1);
%      a1  :  5.786843736434342     -3.341035788988665      0.000000000000000
%      a2  :  0.000000000000000      6.682071577977333      0.000000000000000
%      a3  :  0.000000000000000      0.000000000000000      28.467668849960234

%                Wyckoff positions
% Number of Wyckoff positions :    3
% No.   Element    X                      Y                      Z                      Concentration
%    1    Ti       0.000000000000000      0.000000000000000      0.101711091634490      1.000000000000000
%    2    Se       0.333333333333333     -0.333333333333333      0.203422183268980      1.000000000000000
%    3    Se      -0.333333333333333      0.333333333333333      0.000000000000000      1.000000000000000
str2='               Wyckoff positions'
wyck_pos=strncmp(str2,splitLines,length(str2));
wyck_pos_index=find(wyck_pos==1)
%wyck_pos_index=484;%input('line '    );
%wyck_pos_index=553;
wyck_pos_index=596
element=true;
inde=0;
while element
    inde=inde+1;
wyck_str=splitLines{wyck_pos_index+inde};
         pos_element=strfind(wyck_str,num2str(inde));
         element=~isempty(pos_element);
         if element
         data=str2num(wyck_str(13:end));
         wyck{inde}.pos=data(2:4);
         wyck{inde}.Name='Ti'
         end;
end

      for ind=1:3
         latt_str=splitLines{lat_pos_index1+ind};
         pos=strfind(latt_str,':');
         latt_vec(ind,:)=str2num(latt_str(pos+1:end));
      end;
lat_pos_index=find(lat_pos==1);
abc_str=splitLines{lat_pos_index};
abc=str2num(abc_str(length(str0)+1:end));
    % find number of voxels in each dimension
    str0='END   processing Wannier bands';
    lc_gridpositions=strncmp(str0,splitLines,length(str0));
lc_gridindex=find(lc_gridpositions==1);
lc_gridindex(end+1)=numel(splitLines);
if number== -1
    tbstart=lc_gridindex(end-1);
    tbend=lc_gridindex(end);
else
    tbstart=lc_gridindex(number);
    tbend=lc_gridindex(number+1);
end;
lines=splitLines(tbstart:tbend);
str0='WF(';
cl={str0}
lc_gridpositions=strncmp(str0,lines,length(str0));
orb1_orb2_index=cellfun(@strfind,lines,repmat(cl,numel(lines),1),'UniformOutput',false);
startindices=find(cellfun(@isempty,orb1_orb2_index)==0);
str0='------------- Grid definition ------------';
endpos=strncmp(str0,lines,length(str0));
endpos_index=find(endpos==1);
startindices(end+1)=endpos_index(1)-2;
wannierstrings={};
tb=zeros(numel(lines),6+complx);
    pos=0;
for num=1:numel(startindices)-1
    % iterate over all combinations of wannier functions
    % get name of WF1 and WF2
    str0=lines{startindices(num)};
    startlabel=regexp(str0,'\(');
    endlabel=regexp(str0,'\)');
    label1=str0(startlabel(1)+1:endlabel(1)-1)
    label2=str0(startlabel(2)+1:endlabel(2)-1)
    % find the labels in a set of labels
     orb1_f=strncmp(label1,wannierstrings,length(str0));
     orb1=find(orb1_f==1);
     if isempty(orb1)
         wannierstrings{end+1}=label1;
         orb1=numel(wannierstrings);
     end;
     orb2_f=strncmp(label2,wannierstrings,length(str0));
     orb2=find(orb2_f==1);
     if isempty(orb2)
         wannierstrings{end+1}=label2;
         orb2=numel(wannierstrings);
     end;
     % fill in to tb
     startlabel=cellfun(@regexp,lines(startindices(num)+1:startindices(num+1)-1),repmat({'='},startindices(num+1)-1-(startindices(num)),1),'UniformOutput' ,false);
     for s=1:numel(startlabel)
         strng=lines{startindices(num)+s}
         try
             r=str2num(strng(startlabel{s}(1)+1:startlabel{s}(2)-4));
         catch
             r=''
         end
         try
         t=str2num(strng(startlabel{s}(2)+1:end));
         t=real(t)+imag(t);
         pos=pos+1;
         tb(pos,:)=[r,orb1,orb2,t(1:(1+complx))];
         catch
             t
         end
         
     end
end;
tb=tb(1:pos,:);
if false
tb(:,1)=round(tb(:,1)/abc(1)*10000)/10000;
tb(:,2)=round(tb(:,2)/abc(2)*10000)/10000;
tb(:,3)=round(tb(:,3)/abc(3)*10000)/10000;
else
% some post processing (does not work automatically, need to set the ranges
% etc correctly for the model
if false
    wyck{1}.range=1:5;
    wyck{2}.range=6:8;
    wyck{3}.range=9:11;
    orbs=[];
    for ned=1:3
        orbs=[orbs,ned*ones(1,size(wyck{ned}.range,2))];
    end;
    for n=1:size(tb,1)
        atom1=orbs(tb(n,4));
        atom2=orbs(tb(n,5));
        tb(n,1:3)=round((tb(n,1:3)+(wyck{atom1}.pos-wyck{atom2}.pos))*10000)/10000;
    end;
end
% do a transformation to the basis of the lattice vectors
tb1=zeros(size(tb,1),3);
for n=1:size(tb,1)
   tb1(n,:)= round(tb(n,1:3)/latt_vec*10000)/10000;
end;
end


tb(:,1:3)= tb1;
dlmwrite([outfile,'model_',num2str(number),'.csv'],tb,'delimiter',',','precision',7);
toc

end