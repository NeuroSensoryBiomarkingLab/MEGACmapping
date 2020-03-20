clear sound;

%Loads the noise that will be reproduced the whole run
noise = dir('*Noise_run*');
noise = noise.name;
runNumb = noise(10);
load(noise);
x = Gated_noiseVector(:);
    sound(x,44100);

%Loads the audio to be played
stim = dir('*SparsedRandTonesStimuli_run*');
stim = stim.name;
load(stim);

[rows, cols] = size(TOTAL_MAT_file_reduced);

for cell = 1:rows
    
    WaitingTime = KbTriggerWait(53);
    WaitingTime_List{cell} = [cell WaitingTime]; %#ok<SAGROW>
    
    y = TOTAL_MAT_file_reduced(cell,:);
    sound(y,44100);
    
    %%function secs = KbTriggerWait(keyCode, deviceNumber)
    
    %%Add a silence. If the fMRI acquires just after the trigger, the
    %%stimuli should be done with the silence before the tones.

end

WaitingTime_File = cell2mat(WaitingTime_List);




dlmwrite(['SparsedRandTonesStimuli_WaitingTimeList_run' num2str(runNumb) '.txt'],WaitingTime_File,'delimiter','\t', 'precision', 9);

fprintf('The run is finished.\n');

clear sound;
