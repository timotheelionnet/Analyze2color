destDir = 'code';
sourceDir = 'oldCode';
addpath('oldCode');
fList = matlab.codetools.requiredFilesAndProducts([sourceDir,'/Analyze2color_diatrack.m']);
fList = [fList, matlab.codetools.requiredFilesAndProducts([sourceDir,'/Analyze2color_diatrack2.m'])];
fList = [fList, matlab.codetools.requiredFilesAndProducts([sourceDir,'/Analyze2color_diatrack3.m'])];
fList = [fList, matlab.codetools.requiredFilesAndProducts([sourceDir,'/Analyze2color_diatrack4.m'])];
fList = [fList, matlab.codetools.requiredFilesAndProducts([sourceDir,'/Analyze2color_diatrack5.m'])];

fList = unique(fList)';

for i=1:numel(fList)
    
    destPath = strrep(fList{i},sourceDir,destDir);
    copyfile(fList{i},destPath);
end