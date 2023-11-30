function f=ConvPosFileTxtToMat(posFolderName)
% posFolderName should be the folder containing only position text files
cd(posFolderName);
fullfilenames=dir('*.txt');
for i=1:length(fullfilenames)
    [filepath,filename,fileext]=fileparts(fullfilenames(i).name);
    r=load([filename,'.txt']);
    save([filename,'.mat'],'r');
end
end