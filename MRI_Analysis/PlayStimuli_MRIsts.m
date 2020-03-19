clear sound;

%Loads the noise that will be reproduced the whole run
load Noise_run4.mat;
x = Gated_noiseVector(:);
    sound(x,44100);

%Loads the audio to be played
load SparsedRandTonesStimuli_run4.mat;

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

dlmwrite('SparsedRandTonesStimuli_run4_WaitingTimeList.txt',WaitingTime_File,'delimiter','\t', 'precision', 9);

fprintf('The run is finished.\n');

clear sound;