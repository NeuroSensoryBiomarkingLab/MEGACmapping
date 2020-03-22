clear all

if isempty(dir('*freq0.txt'))
    fprintf('freq0.txt not present. I will import .txt files.\n')
else
    mkdir freq0
    movefile('*freq0.txt', 'freq0')
    fprintf('As freq0.txt file was present, I moved it to folder <freq0>.\n')
end

currentDir = dir('Run*.txt');
fileNames = {currentDir(:,1).name}';
onsets={};
names={};
durations={};
Number_Run_files=length(count(fileNames,"Run"));
ShortFileNames={};


%% Get list of files for Onset
if Number_Run_files<1
    error("There isn't any Run file in this folder");
else
    for jj=1:length(fileNames)
        
        ShortFileNames{jj}=fileNames{jj}(1:4);
    end
    
    ShortFileNames=unique(ShortFileNames);
    Number_Diff_runs=length(ShortFileNames);
 
    
    for j=1:length(ShortFileNames)
        Run_dir=dir(cell2mat(strcat(ShortFileNames(j),'*.txt')));
        Run_filenames={Run_dir(:,1).name};
        
        for i = 1:length(Run_filenames)
            var=readtable(Run_filenames{i});
            var=table2array(var);
            var1=(var(1:end,1))';
            onsets{i}=var1;
            var2=(var(1:end,2))';
            durations{i}=var2;
            var3=extractBetween(Run_filenames{i},'_','.txt');
            names{i} = cell2mat(var3);
        end
        Run_name = strcat(ShortFileNames(j),'.mat');
        save(char(Run_name), 'names', 'onsets', 'durations');
    end
end