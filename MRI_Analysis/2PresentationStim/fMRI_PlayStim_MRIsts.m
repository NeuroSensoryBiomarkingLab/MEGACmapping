clear sound;

%Loads the noise that will be reproduced the whole run
noise = dir('*Noise_run*,mat');
noise = noise.name;
runNumb = noise(10);
load(noise);
x = Gated_noiseVector(:);
    sound(x,44100);

%Loads the audio to be played
stim = dir('*SparsedRandTonesStimuli_run*.mat');
stim = stim.name;
load(stim);

[rows, cols] = size(TOTAL_MAT_file_reduced);

for cell = 1:rows
    
    WaitingTime = KbTriggerWait(53); %% Waits for the trigger given by the fMRI scanner in the form of a number 5 keyboard press.
    WaitingTime_List{cell} = [cell WaitingTime]; %#ok<SAGROW>
    
    y = TOTAL_MAT_file_reduced(cell,:);
    sound(y,44100);
    
end

WaitingTime_File = cell2mat(WaitingTime_List); %%Creates a txt file with the time between each trigger.

dlmwrite(['SparsedRandTonesStimuli_WaitingTimeList_run' num2str(runNumb) '.txt'],WaitingTime_File,'delimiter','\t', 'precision', 9);

fprintf('The run is finished.\n');

clear sound; %% Clears the noise that is longer than the stim.
