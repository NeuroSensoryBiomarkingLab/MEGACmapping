function varargout = process_STRF_mapping( varargin )

% PROCESS_STRF_MAPPING:
% Input: In 'Process 2', Scout time series (Files A), Raw recording (Files B)
% Output: STRFs using reverse correlation and cortical maps of STRF
% characteristics

% Compatible with Brainstorm version 20-Sep-2017 which is available for download at:
% https://neuroimage.usc.edu/brainstorm/

% To use this script, copy paste it into the following folder:
% ~/.brainstorm/process

% You will find the process window under the category "STRF"

macro_methodcall;
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
% Description the process
sProcess.Comment     = 'STRF Mapping (A = Scout Time Series), B = Raw Data)';
sProcess.FileTag     = '';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'STRF';
sProcess.Index       = 400;

% Label
sProcess.options.label_about.Comment = '<BR><I>Depending on your RAM availability, you may have to divide your scout-time series (Files A) into smaller trials (we recommend trials of 10 seconds)</I>';
sProcess.options.label_about.Type    = 'label';

% Definition of the input accepted by this process
sProcess.InputTypes  = {'raw', 'matrix'};
sProcess.OutputTypes = {'raw', 'data', 'matrix', 'results'};
sProcess.nInputs     = 2;
sProcess.nMinFiles   = 1;
sProcess.isPaired    = 0;
sProcess.isSeparator = 0;

% Channel name
sProcess.options.channelname.Comment = 'Channel name: ';
sProcess.options.channelname.Type    = 'hidden';
sProcess.options.channelname.Value   = 'UADC001';

% Label
sProcess.options.label_input.Comment = '<BR><B><U>Input Options</U></B>:';
sProcess.options.label_input.Type    = 'label';

% Event Name Stem
sProcess.options.eventStem.Comment = 'Stem for tone-pip event names (e.g. if you named events freq_1, freq_2, etc., then stem is freq_): ';
sProcess.options.eventStem.Type    = 'text';
sProcess.options.eventStem.Value   = 'freq_';

% Name for all events group
sProcess.options.allEventsName.Comment = 'Name given to the event group containing all events: ';
sProcess.options.allEventsName.Type    = 'text';
sProcess.options.allEventsName.Value   = 'allTTL';

% Label
sProcess.options.label_description.Comment = '<BR><B><U>Output Options</U></B>:';
sProcess.options.label_description.Type    = 'label';

% Save directory
sProcess.options.saveDir.Comment = 'Save directory for STRFs: ';
sProcess.options.saveDir.Type    = 'text';
sProcess.options.saveDir.Value   = 'X:/----/----';

% Description
sProcess.options.description.Comment = 'Output description (short descriptive text that will form the filename stem of output files): ';
sProcess.options.description.Type    = 'text';
sProcess.options.description.Value   = '';

% Label
sProcess.options.label_zscore.Comment = '<BR><B><U>Spike Detection Options</U></B>:';
sProcess.options.label_zscore.Type    = 'label';

% Z-score threshold
sProcess.options.threshold.Comment = 'Activation z-score threshold:';
sProcess.options.threshold.Type = 'value';
sProcess.options.threshold.Value = {1.0, '', 3};

% Z-score pre-buff
sProcess.options.preBuff.Comment = 'Minimum duration of z-score baseline segments:';
sProcess.options.preBuff.Type = 'value';
sProcess.options.preBuff.Value = {0.100, 'ms', 0};

% Z-score post-buff
sProcess.options.postBuff.Comment = 'Minimum delay after a tone-pip before retaking a new z-score baseline:';
sProcess.options.postBuff.Type = 'value';
sProcess.options.postBuff.Value = {0.350, 'ms', 0};

% Label
sProcess.options.label_strf.Comment = '<BR><B><U>STRF Options</U></B>:';
sProcess.options.label_strf.Type    = 'label';

% STRF length
sProcess.options.strfLength.Comment = 'STRF duration (MULTIPLES OF 100):';
sProcess.options.strfLength.Type = 'value';
sProcess.options.strfLength.Value = {0.500, 'ms', 0};

% STRF bin size
sProcess.options.binSize.Comment = 'STRF bin size:';
sProcess.options.binSize.Type = 'value';
sProcess.options.binSize.Value = {0.004, 'ms', 0};

% STRF smoothing options
sProcess.options.isMedfilt.Comment = 'Use median filtering? (1 = yes, 0 = no)';
sProcess.options.isMedfilt.Type = 'value';
sProcess.options.isMedfilt.Value = {0, '', 0};

sProcess.options.isSmooth.Comment = 'Use gaussian smoothing? (1 = yes, 0 = no)';
sProcess.options.isSmooth.Type = 'value';
sProcess.options.isSmooth.Value = {1, '', 0};

% Baseline duration for calculation of STRF Z-score
sProcess.options.strfBaseline.Comment = 'Duration of baseline at the beginning of STRFs for calculation of z-score:';
sProcess.options.strfBaseline.Type = 'value';
sProcess.options.strfBaseline.Value = {0.150, 'ms', 0};

% Reference for calculation of bandwidth from peak
sProcess.options.baseDist.Comment = 'Standard deviation from gaussian-STRF peak to measure bandwidth (FWHM ~ 2.355):';
sProcess.options.baseDist.Type = 'value';
sProcess.options.baseDist.Value = {2.355, '', 3};

% M50 Window
sProcess.options.rM50Window.Comment = 'M50-Peak Search Window:';
sProcess.options.rM50Window.Type = 'range';
sProcess.options.rM50Window.Value = {[0.030, 0.075], 'ms', 0};

% M100 Window
sProcess.options.rM100Window.Comment = 'M100-Peak Search Window:';
sProcess.options.rM100Window.Type = 'range';
sProcess.options.rM100Window.Value = {[0.075, 0.135], 'ms', 0};

% M50 Window for Gaussian
sProcess.options.gM50Window.Comment = 'M50 Gaussian Fit Window:';
sProcess.options.gM50Window.Type = 'range';
sProcess.options.gM50Window.Value = {[0.020, 0.085], 'ms', 0};

% M100 Window for Gaussian
sProcess.options.gM100Window.Comment = 'M100 Gaussian Fit Window:';
sProcess.options.gM100Window.Type = 'range';
sProcess.options.gM100Window.Value = {[0.075, 0.150], 'ms', 0};
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
Comment = 'STRF Mapping';
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputsA, sInputsB)

% Setup progress bar
progressPos = bst_progress('get');
bst_progress('text', 'Setting up variables...');
    
% Setup process window options
OPTIONS.eventStem = strtrim(sProcess.options.eventStem.Value);
OPTIONS.allEventsName = strtrim(sProcess.options.allEventsName.Value);
OPTIONS.saveDir = strtrim(sProcess.options.saveDir.Value);
OPTIONS.description = strtrim(sProcess.options.description.Value);
OPTIONS.threshold  = sProcess.options.threshold.Value{1};
OPTIONS.preBuff = sProcess.options.preBuff.Value{1}*1000;
OPTIONS.postBuff = sProcess.options.postBuff.Value{1}*1000;
OPTIONS.strfLength = sProcess.options.strfLength.Value{1}*1000;
OPTIONS.binSize = sProcess.options.binSize.Value{1}*1000;
OPTIONS.isMedfilt = sProcess.options.isMedfilt.Value{1};
OPTIONS.isSmooth = sProcess.options.isSmooth.Value{1};
OPTIONS.strfBaseline = sProcess.options.strfBaseline.Value{1}*1000;
OPTIONS.baseDist = sProcess.options.baseDist.Value{1};
OPTIONS.rM50Window = sProcess.options.rM50Window.Value{1}*1000;
OPTIONS.rM100Window = sProcess.options.rM100Window.Value{1}*1000;
OPTIONS.gM50Window = sProcess.options.gM50Window.Value{1}*1000;
OPTIONS.gM100Window = sProcess.options.gM100Window.Value{1}*1000;

% Setup save directory for output STRF files and databases
mkdir(strcat(OPTIONS.saveDir, '/', sInputsB.SubjectName), sInputsA(1).Condition);
mkdir(strcat(OPTIONS.saveDir, '/', sInputsB.SubjectName, '/', sInputsA(1).Condition), OPTIONS.description);

% Get raw file data
rawMat = in_bst_data(sInputsB(1).FileName);
% Get all events
events = rawMat.F.events;
iFreqEvents = find(contains({events.label}, OPTIONS.eventStem));
totalTTL(1:length(iFreqEvents)) = deal(0);
% Get list of bad segments in raw file
badSeg = panel_record('GetBadSegments', rawMat.F);
badSeg = badSeg+1; %Brainstorm indexing starts at 0, so add 1 to make it start at 1
% Get ChannelFlag
ChannelFlag = rawMat.ChannelFlag;

% Load channel file
ChannelMat = in_bst_channel(sInputsB(1).ChannelFile);
% Get channel to process
chanName = strtrim(sProcess.options.channelname.Value);
iChannel = find(strcmpi({ChannelMat.Channel.Name}, chanName));
% Read channel to process
SamplesBounds = [];
[~, TimeVector] = in_fread(rawMat.F, ChannelMat, 1, SamplesBounds, iChannel);

% For each file A (scout time series), detect significant activations with
% detectSpikes
for iFile = 1:length(sInputsA)
    
    % Progress bar
    bst_progress('text', 'Detecting spikes...');
    
    % Get scout time series
    scoutMat = in_bst_data(sInputsA(iFile).FileName);
    
    if iFile == 1
        % Get list of scout vertices
        scoutVertices = scoutMat.Atlas.Scouts.Vertices;
        % Initialise threshVals
        threshVals(1:length(scoutVertices)) = struct('Values', [], 'TimeValues', []);
        % Initialize empty variables
        lastData = [];
        lastMask = [];
        lastBaseline = [];
    end
    
    % Detect significant activations using detectSpikes
    [vals, lastD, lastM, lastB, totTTL] = detectSpikes(scoutMat, badSeg, OPTIONS, lastData, lastMask, lastBaseline, totalTTL, threshVals, TimeVector, events, iFreqEvents);
    threshVals = vals;
    lastData = lastD;
    lastMask = lastM;
    lastBaseline = lastB;
    totalTTL = totTTL;
    
    % Update progress bar
    bst_progress('set', progressPos + round((iFile / length(sInputsA)) * 100));
end

% Save threshVals which represents the threshold values for each
% significant activation
bst_progress('text', 'Saving activation values...');
save(strcat(OPTIONS.saveDir, '/', sInputsB.SubjectName, '/', sInputsA(1).Condition, '/', OPTIONS.description, '/threshVals.mat'), 'threshVals');

% Get protocol info
ProtocolInfo = bst_get('ProtocolInfo');

% Generate maps
mapGen(ChannelFlag, threshVals, scoutVertices, totalTTL, OPTIONS, sInputsA(1).iStudy, sInputsA(1).Condition, sInputsB.SubjectName, scoutMat.SurfaceFile, ProtocolInfo, events, iFreqEvents);

% Return the input files
OutputFiles = {sInputsB(1).FileName};
end


%% ===== Spike Detection =====
function [vals, lastD, lastM, lastB, totTTL] = detectSpikes(scoutMat, badSeg, OPTIONS, lastData, lastMask, lastBaseline, totalTTL, threshVals, TimeVector, events, iFreqEvents)

% Create Bad Segment Mask
toRemove = false(1,size(badSeg,2));
if ~isempty(badSeg)
    for ii = 1:size(badSeg,2)
        startBad = TimeVector(badSeg(1,ii));
        endBad = TimeVector(badSeg(2,ii));
        if ~isempty(find(scoutMat.Time == startBad, 1)) && ~isempty(find(scoutMat.Time == endBad, 1)) %If the whole bad segment is contained within this trial
            badSeg(1,ii) = find(scoutMat.Time == startBad);
            badSeg(2,ii) = find(scoutMat.Time == endBad);
        elseif isempty(find(scoutMat.Time == startBad, 1)) && ~isempty(find(scoutMat.Time == endBad, 1)) %If the start of the bad segment begins before this trial
            badSeg(1,ii) = 1;
            badSeg(2,ii) = find(scoutMat.Time == endBad);
        elseif ~isempty(find(scoutMat.Time == startBad, 1)) && isempty(find(scoutMat.Time == endBad, 1)) %If the end of the bad segment ends after this trial
            badSeg(1,ii) = find(scoutMat.Time == startBad);
            badSeg(2,ii) = length(scoutMat.Time);
        else % If bad segment is not within this trial
            toRemove(ii) = true;
        end
    end
    badSeg(:,toRemove) = [];
    % Create scout time series mask
    Smask = true(size(scoutMat.Time));
    % Loop on each segment: mark as bad
    for iSeg = 1:size(badSeg, 2)
        Smask(badSeg(1,iSeg):badSeg(2,iSeg)) = false;
    end
end

% Calculate sampling frequency
sFreq = round( 1 ./ (scoutMat.Time(2) - scoutMat.Time(1)));

% Extract the timeSeries and take the absolute value...should already be
% meancorrected from DC offset removal
timeSeries = abs(scoutMat.Value);

%Find all TTL event times in this scout time series
iallTTL = find(strcmp({events.label}, OPTIONS.allEventsName));
iTTLevents = find(events(iallTTL).times >= scoutMat.Time(1) & events(iallTTL).times <= scoutMat.Time(end));

%Convert the time value to the nearest sample
TTLSamples = zeros(1, length(iTTLevents));
for iii = 1:length(iTTLevents)
    [~, iMin]= min(abs(scoutMat.Time - events(iallTTL).times(iTTLevents(iii)))); %find difference between target and vector and then the minimum to find nearest value
    TTLSamples(iii) = iMin; %this is the nearest timepoint to the event
end

% Find total number of tone pips for each frequency within this scout time
% series and add this number to totalTTL to get a total count of TTLs
% analysed for all input files A
for iFreq = 1:length(iFreqEvents)
    iTrialFreqEvents = find(events(iFreqEvents(iFreq)).times >= scoutMat.Time(1) & events(iFreqEvents(iFreq)).times <= scoutMat.Time(end));
    
    %Remove events that fall into bad segments
    toRemove = false(1,length(iTrialFreqEvents));
    for ii = 1:length(iTrialFreqEvents)
        if Smask(scoutMat.Time == events(iFreqEvents(iFreq)).times(iTrialFreqEvents(ii)))
            continue;
        else
            toRemove(ii) = true;
        end
    end
    iTrialFreqEvents(toRemove)=[];
    
    if ~isempty(iTrialFreqEvents)
        totalTTL(iFreq) = totalTTL(iFreq) + length(iTrialFreqEvents);
    end
end
totTTL = totalTTL;

% DYNAMIC ZSCORE TRANSFORMATION OF SCOUT TIME SERIES
if ~isempty(iTTLevents)
    % Set buffer durations
    preBuff = round(-1*(OPTIONS.preBuff/1000)*sFreq); %in milliseconds, then convert into samples
    postBuff = round((OPTIONS.postBuff/1000)*sFreq); %in milliseconds, then convert into samples
    
    % Loop through each TTL event and determine if a new z-score baseline
    % is calculated
    for ii = 1:length(TTLSamples) 
        if ii == 1 && ~isempty(lastBaseline) % Start by calculating the zscore for everything up to the first event using the last baseline
            initSegData = timeSeries(:, 1:round(TTLSamples(ii)-1));
            initBaselineMean = mean(lastBaseline, 2);
            initBaselineStd = std(lastBaseline, 0, 2);
            
            % Transform the mean and std into a matrix the same size as the
            % segData in order to perform the subtraction
            initBaselineMean = repmat(initBaselineMean, 1, size(initSegData, 2));
            initBaselineStd = repmat(initBaselineStd, 1, size(initSegData, 2));
            initSegZscore = (initSegData - initBaselineMean)./initBaselineStd;
            zscoreTimeSeries = initSegZscore;
        end
        
        if ii == 1 && isempty(lastBaseline) % If this is the first block and the first event, take whatever is the initial segment baseline
            baselineData = timeSeries(:, 1:TTLSamples(ii)-1);
            
            % Check if part of the baseline includes a BAD segment and remove that part
            isBad = find(Smask(1:TTLSamples(ii)-1) == 0);
            if ~isempty(isBad)
                baselineData(:,isBad) = [];
            end
            zscoreTimeSeries = zeros(size(timeSeries,1), TTLSamples(ii)-1); %Add a filler at the front to keep same size as timeSeries
            
        else % For everything after the initSeg of the scoutTimeSeries
            if ii == 1 % If this is the first event, use the lastData as a starting point to determine if you
                if TTLSamples(ii)-1 + size(lastData,2) >= (postBuff + abs(preBuff)) % If there is enough space, retake a new baseline
                    if size(lastData, 2) > postBuff % If part of the baseline can start in lastData, use it
                        baselineData = cat(2, lastData(:, round(postBuff+1):end), timeSeries(:, 1:round(TTLSamples(ii)-1)));
                        
                        % Check if part of the baseline includes a BAD segment and remove that part
                        isBad = find(lastMask(round(postBuff+1):end) == 0);
                        isBad = cat(2, isBad, (find(Smask(1:round(TTLSamples(ii)-1)) == 0)) + round(length(lastMask)-postBuff)); %Readjust indexes so that the follow the lastMask ones
                        if ~isempty(isBad)
                            baselineData(:,isBad) = [];
                        end
                        
                    else % Only use the current scoutTimeSeries, starting after when lastData's postBuff would have ended
                        baselineData = timeSeries(:, round(postBuff - size(lastData,2) + 1):TTLSamples(ii)-1);
                        
                        % Check if part of the baseline includes a BAD segment and remove that part
                        isBad = find(Smask(round(postBuff - size(lastData,2) + 1):TTLSamples(ii)-1) == 0);
                        if ~isempty(isBad)
                            baselineData(:,isBad) = [];
                        end
                    end
                    
                else  % If not enough space for new baseline even with lastData, use previous baseline
                    baselineData = lastBaseline;
                end
                
            else % For all events after the first event
                if TTLSamples(ii) - TTLSamples(ii-1) >= (postBuff + abs(preBuff)) %Only if theres more than prebuff+post-buff between events do you reset the baselineData…if not it will just use the previous baselineData.
                    baselineData = timeSeries(:, round(TTLSamples(ii-1) + postBuff): round(TTLSamples(ii)-1));
                    
                    % Check if part of the baseline includes a BAD segment and remove that part
                    isBad = find(Smask(round(TTLSamples(ii-1) + postBuff): round(TTLSamples(ii)-1)) == 0);
                    if ~isempty(isBad)
                        baselineData(:,isBad) = [];
                    end
                end
            end
        end
        
        if isempty(baselineData) && ii == 1 % If there is no baselineData that lies outside of a BAD segment, use the last one
            if ~isempty(lastBaseline)
                baselineData = lastBaseline;
            else
                segData = timeSeries(:, TTLSamples(ii):round(TTLSamples(ii+1)-1)); % Take the segment from the current event to the sample just before the next event
                zscoreTimeSeries = cat(2, zscoreTimeSeries, zeros(size(segData,1),size(segData,2))); % Add zeros until the next TTL event
                continue;
            end
        elseif isempty(baselineData) % In the case where this is no longer the first TTL
            if isempty(previousBaseline)
                segData = timeSeries(:, TTLSamples(ii):round(TTLSamples(ii+1)-1)); % Take segment from the current event to the sample just before the next event
                zscoreTimeSeries = cat(2, zscoreTimeSeries, zeros(size(segData,1),size(segData,2))); % Add zeros until the next TTL event
                continue;
            else
                baselineData = previousBaseline;
            end
        end
        
        if ii == length(TTLSamples) % If this is the last event
            lastD = timeSeries(:, TTLSamples(ii):end);
            lastM = Smask(TTLSamples(ii):end);
            lastB = baselineData;
            segData = timeSeries(:,TTLSamples(ii):end);
        else
            segData = timeSeries(:, TTLSamples(ii):round(TTLSamples(ii+1)-1)); % Take segment from the current event to the sample just before the next event
        end
        
        baselineMean = mean(baselineData, 2);
        baselineStd = std(baselineData, 0, 2);
        
        % Transform the mean and std into a matrix the same size as the
        % segData in order to perform the subtraction
        baselineMean = repmat(baselineMean, 1, size(segData, 2));
        baselineStd = repmat(baselineStd, 1, size(segData, 2));
        
        % Transform time series into z-score
        segZscore = (segData - baselineMean)./baselineStd;
        zscoreTimeSeries = cat(2, zscoreTimeSeries, segZscore);
        previousBaseline = baselineData;
    end

    % Remove BAD segments from z-score time series so that no spikes can be
    % detected in BAD segments
    zscoreTimeSeries(:,~Smask) = deal(0);
    
else % If there is no TTL in the segment, restart from scratch is if it was the first segment
    lastD = [];
    lastM = [];
    lastB = [];
    zscoreTimeSeries = zeros(size(timeSeries,1), size(timeSeries,2));
end

% FIND SIGNIFICANT Z-SCORE ACTIVATIONS
for iVertex = 1:size(zscoreTimeSeries,1)
    
    % First find local maxima
    [actVals, iActVals] = findpeaks(zscoreTimeSeries(iVertex,:));
    
    % Save only the peaks that cross threshold
    iActKeep = find(actVals > OPTIONS.threshold);
    actVals = actVals(iActKeep);
    iActVals = iActVals(iActKeep);
    
    % Concatenate actVals and iActVals into one matrix
    if ~isempty(iActVals) %&& ~isempty(iInhibVals)
        allVals(1,:) = iActVals; % First row = index where an event occurs
        allVals(2,:) = actVals; % Second row = value of the event
    else % If both are empty, move on to next iVertex
        allVals = [];
        continue;
    end
    
    % Sort actVals
    allVals = transpose(allVals);
    allVals = sortrows(allVals,1); % Sort according to indices (stored in column 1)
    allVals = transpose(allVals); % Transpose back into row form
    threshVals(iVertex).Values = cat(2, threshVals(iVertex).Values, allVals(2,:));
    threshVals(iVertex).TimeValues = cat(2, threshVals(iVertex).TimeValues, scoutMat.Time(allVals(1,:)));
    allVals = [];
end

vals = threshVals;

end


%% ===== MAP GENERATION =====
function [] = mapGen(ChannelFlag, threshVals, vertices, totalTTL, OPTIONS, iStudy, conditionName, subjectName, surfaceFileName, ProtocolInfo, events, iFreqEvents)

% Set progress bar
bst_progress('text', 'Generating STRFs and maps...');
bst_progress('set', 0);

%% PREPARE OUTPUT FILES
sourceStem = OPTIONS.description;
figureFileStem = strcat(OPTIONS.saveDir, '/', subjectName, '/', conditionName, '/', OPTIONS.description, '/');
mapFileStem = strcat(ProtocolInfo.STUDIES, '/', subjectName, '/', conditionName, '/results_STRF_mapping_', char(datetime('now','TimeZone','local','Format','yyMMddHHmm')) , '_');

% Initialize source files for STRF maps
iMap = 1;
[M100_filename_peak, iMap] = filename_gen(mapFileStem, iMap);
[M50_filename_peak, iMap] = filename_gen(mapFileStem, iMap);
[M100_filename_zscore, iMap] = filename_gen(mapFileStem, iMap);
[M50_filename_zscore, iMap] = filename_gen(mapFileStem, iMap);
[M100_filename_lat, iMap] = filename_gen(mapFileStem, iMap);
[M50_filename_lat, iMap] = filename_gen(mapFileStem, iMap);

% Initialize source files for gaussian-fitted STRF maps
[M100_filename_peak_GAUSSIAN, iMap] = filename_gen(mapFileStem, iMap);
[M50_filename_peak_GAUSSIAN, iMap] = filename_gen(mapFileStem, iMap);
[M100_filename_band_GAUSSIAN, iMap] = filename_gen(mapFileStem, iMap);
[M50_filename_band_GAUSSIAN, iMap] = filename_gen(mapFileStem, iMap);
[M100_filename_lat_GAUSSIAN, iMap] = filename_gen(mapFileStem, iMap);
[M50_filename_lat_GAUSSIAN, iMap] = filename_gen(mapFileStem, iMap);
[M100_filename_tempmod_GAUSSIAN, iMap] = filename_gen(mapFileStem, iMap);
[M50_filename_tempmod_GAUSSIAN, ~] = filename_gen(mapFileStem, iMap);

% Initialize STRF and gaussian-fitted STRF feature databases
STRF_features.columns = {'vertex', 'peak', 'lat', 'zscore'};
STRF_features.M50 = zeros(length(threshVals), length(STRF_features.columns));
STRF_features.M100 = zeros(length(threshVals), length(STRF_features.columns));
GAUSSIAN_features.columns = {'vertex', 'peak', 'bandwidth', 'latency', 'temporal modulation'};
GAUSSIAN_features.M50 = zeros(length(threshVals), length(GAUSSIAN_features.columns));
GAUSSIAN_features.M100 = zeros(length(threshVals), length(GAUSSIAN_features.columns));

% Prep the STRF/Gaussian matrix database
STRF = zeros(length(iFreqEvents), round(OPTIONS.strfLength/OPTIONS.binSize), length(threshVals));
M50_gaussian = zeros(length(iFreqEvents), round(OPTIONS.strfLength/OPTIONS.binSize), length(threshVals));
M100_gaussian = zeros(length(iFreqEvents), round(OPTIONS.strfLength/OPTIONS.binSize), length(threshVals));

%Get event frequency values and save them for reference
freqVals = cell(1,length(iFreqEvents));
for ii = 1:length(iFreqEvents)
    freqString = fliplr(events(iFreqEvents(ii)).label); %Flip the string so the vertex number appears first
    freqString = fliplr(strtok(freqString, '_')); % flip it back once you have just the number before the _
    freqVals{ii} = freqString;
end
save(strcat(figureFileStem, 'freqVals.mat'), 'freqVals');

% Save scout vertices
save(strcat(figureFileStem, 'allVertices.mat'), 'vertices');

% Find how many vertices in the cortical surface
cortexVertices = load(strcat(ProtocolInfo.SUBJECTS, '/', surfaceFileName), 'Vertices');
totalVertices = length(cortexVertices.Vertices);

% Get head model filename
hm = bst_get('HeadModelForStudy', iStudy);
hm_filename = hm.FileName;
hm = [];

% Create source structure
source.ImagingKernel = [];
source.ImageGridAmp = zeros(totalVertices, 1);
source.Std = [];
source.Whitener = [];
source.nComponents = 1;
source.Comment = OPTIONS.description;
source.Function = 'wmne';
source.Time = 1;
source.DataFile = '';
source.HeadModelFile = hm_filename;
source.HeadModelType = 'surface';
source.ChannelFlag =  ChannelFlag;
source.GoodChannel = [];
source.SurfaceFile = surfaceFileName;
source.Atlas = [];
source.GridLoc = [];
source.GridOrient = [];
source.GridAtlas = [];
source.Options = [];
source.ColormapType = [];
source.ZScore = [];
source.nAvg = 1;
source.displayUnits = 1;
source.History = [];

% Initialize all sources used for mapping STRF features
sourceM50_peak = source;
sourceM100_peak = source;
sourceM50_zscore = source;
sourceM100_zscore = source;
sourceM50_lat = source;
sourceM100_lat = source;
sourceM50_peak.Comment = strcat(sourceStem, '-M50-peak');
sourceM100_peak.Comment = strcat(sourceStem, '-M100-peak');
sourceM50_zscore.Comment = strcat(sourceStem, '-M50-zscore');
sourceM100_zscore.Comment = strcat(sourceStem, '-M100-zscore');
sourceM50_lat.Comment = strcat(sourceStem, '-M50-lat');
sourceM100_lat.Comment = strcat(sourceStem, '-M100-lat');

% Initialize all sources used for mapping gaussian-fitted STRF features
sourceM50_peak_GAUSSIAN = source;
sourceM100_peak_GAUSSIAN = source;
sourceM50_band_GAUSSIAN = source;
sourceM100_band_GAUSSIAN = source;
sourceM50_lat_GAUSSIAN = source;
sourceM100_lat_GAUSSIAN = source;
sourceM50_tempMod_GAUSSIAN = source;
sourceM100_tempMod_GAUSSIAN = source;
sourceM50_peak_GAUSSIAN.Comment = strcat(sourceStem, '-M50-peak-GAUSSIAN');
sourceM100_peak_GAUSSIAN.Comment = strcat(sourceStem, '-M100-peak-GAUSSIAN');
sourceM50_band_GAUSSIAN.Comment = strcat(sourceStem, '-M50-bandSTD-GAUSSIAN');
sourceM100_band_GAUSSIAN.Comment = strcat(sourceStem, '-M100-bandSTD-GAUSSIAN');
sourceM50_lat_GAUSSIAN.Comment = strcat(sourceStem, '-M50-lat-GAUSSIAN');
sourceM100_lat_GAUSSIAN.Comment = strcat(sourceStem, '-M100-lat-GAUSSIAN');
sourceM50_tempMod_GAUSSIAN.Comment = strcat(sourceStem, '-M50-tempMod-GAUSSIAN');
sourceM100_tempMod_GAUSSIAN.Comment = strcat(sourceStem, '-M100-tempMod-GAUSSIAN');

% Calculate compensation (compens1) to correct for slight imbalance in number of TTLs
% per frequency
targetNum = mean(totalTTL);
compens1 = zeros(1,length(totalTTL));
for ii = 1:length(totalTTL)
    compens1(ii) = targetNum/totalTTL(ii);
end

% Initialize counters
has_M50 = false(1, length(threshVals));
has_M100 = false(1, length(threshVals));
has_M50_gaussian = false(1, length(threshVals));
has_M100_gaussian = false(1, length(threshVals));

for iVertex = 1:length(threshVals)
    
    % Generate STRF
    strfMatrix = strfGen(compens1, iFreqEvents, OPTIONS, threshVals(iVertex), events);
    
    % Save STRF matrix
    STRF(:,:,iVertex) = strfMatrix;
   
    %Prep STRF zscore baseline
    baseline = strfMatrix(:,1:round(OPTIONS.strfBaseline/OPTIONS.binSize));
    baseline = baseline(:);
    meanMatrix = mean(baseline);
    stdMatrix = std(baseline,1,1);
    
    % Extract X and Y
    [Ylim, Xlim] = size(strfMatrix);
    Yincrement = length(totalTTL)/Ylim; 
    Y = Yincrement:Yincrement:length(totalTTL);
    Xincrement = OPTIONS.strfLength/Xlim;
    X = -OPTIONS.strfLength:Xincrement:-Xincrement;

    %% STRF Features
    
    % Extract M50 response
    subINT_M50 = strfMatrix(:, round((OPTIONS.strfLength-OPTIONS.rM50Window(2))/Xincrement):round((OPTIONS.strfLength-OPTIONS.rM50Window(1))/Xincrement)); %Only look for peak from -75 to -30ms    
    
    % Extract M50 features
    [has_response, peakloc_col, peakloc_row, zscoreVal] = STRF_feat(subINT_M50, meanMatrix, stdMatrix, OPTIONS.strfLength, OPTIONS.rM50Window(2), Xincrement);
    
    if has_response
        has_M50(iVertex) = true;
        
        % Update Cortex Map
        sourceM50_peak.ImageGridAmp(vertices(iVertex),1) = Y(peakloc_row);
        sourceM50_zscore.ImageGridAmp(vertices(iVertex),1) = zscoreVal;
        sourceM50_lat.ImageGridAmp(vertices(iVertex),1) = X(round(peakloc_col));
        
        %Save values in a seperate matrix
        STRF_features.M50(iVertex,:) = [vertices(iVertex), Y(peakloc_row), X(round(peakloc_col)), zscoreVal];  
    end
    
    
    % Extract M100 response
    subINT_M100 = strfMatrix(:, round((OPTIONS.strfLength-OPTIONS.rM100Window(2))/Xincrement):round((OPTIONS.strfLength-OPTIONS.rM100Window(1))/Xincrement)); %Only look for peak from -75 to -30ms    
    
    % Extract M100 features
    [has_response, peakloc_col, peakloc_row, zscoreVal] = STRF_feat(subINT_M100, meanMatrix, stdMatrix, OPTIONS.strfLength, OPTIONS.rM100Window(2), Xincrement);
    
    if has_response
        has_M100(iVertex) = true;
        
        % Update Cortex Map
        sourceM100_peak.ImageGridAmp(vertices(iVertex),1) = Y(peakloc_row);
        sourceM100_zscore.ImageGridAmp(vertices(iVertex),1) = zscoreVal;
        sourceM100_lat.ImageGridAmp(vertices(iVertex),1) = X(round(peakloc_col));
        
        %Save values in a seperate matrix
        STRF_features.M100(iVertex,:) = [vertices(iVertex), Y(peakloc_row), X(round(peakloc_col)), zscoreVal];  
    end
    

    %% Gaussian-fit STRF Features
    
    % Generate gaussian-fitted STRF
    [M50_g, std_gaussianM50, M100_g, std_gaussianM100] = gaussianGen(strfMatrix, totalTTL, OPTIONS);
    
    % Extract M50 features
    [has_response, peakloc_col, peakloc_row, strfBand, strfTempMod] = GAUSSIAN_feat(M50_g, std_gaussianM50, X, Y, OPTIONS, freqVals);
    
    if has_response
        has_M50_gaussian(iVertex) = true;
        
        % Update Cortex Map
        sourceM50_peak_GAUSSIAN.ImageGridAmp(vertices(iVertex),1) = Y(peakloc_row);
        sourceM50_lat_GAUSSIAN.ImageGridAmp(vertices(iVertex),1) = X(round(peakloc_col));
        sourceM50_band_GAUSSIAN.ImageGridAmp(vertices(iVertex),1) = strfBand;
        sourceM50_tempMod_GAUSSIAN.ImageGridAmp(vertices(iVertex),1) = strfTempMod; 
    
        %Save values in a seperate matrix
        GAUSSIAN_features.M50(iVertex,:) = [vertices(iVertex), Y(peakloc_row), strfBand, X(round(peakloc_col)), strfTempMod];
        
        %Save STRF Gaussian
        M50_gaussian(:,:,iVertex) = M50_g;
    end
    
    
    % Extract M100 features
    [has_response, peakloc_col, peakloc_row, strfBand, strfTempMod] = GAUSSIAN_feat(M100_g, std_gaussianM100, X, Y, OPTIONS, freqVals);
    
    if has_response
        has_M100_gaussian(iVertex) = true;
        
        % Update Cortex Map
        sourceM100_peak_GAUSSIAN.ImageGridAmp(vertices(iVertex),1) = Y(peakloc_row);
        sourceM100_lat_GAUSSIAN.ImageGridAmp(vertices(iVertex),1) = X(round(peakloc_col));
        sourceM100_band_GAUSSIAN.ImageGridAmp(vertices(iVertex),1) = strfBand;
        sourceM100_tempMod_GAUSSIAN.ImageGridAmp(vertices(iVertex),1) = strfTempMod; 
    
        %Save values in a seperate matrix
        GAUSSIAN_features.M100(iVertex,:) = [vertices(iVertex), Y(peakloc_row), strfBand, X(round(peakloc_col)), strfTempMod];
        
        %Save STRF Gaussian
        M100_gaussian(:,:,iVertex) = M100_g;
    end

    bst_progress('set', round((iVertex / length(vertices)) * 100));
end

% Only keep vertex entries in the datamatrices for the vertices where the M50s and M100s were detected
STRF_features.M50(~has_M50,:) = [];
STRF_features.M100(~has_M100,:) = [];
GAUSSIAN_features.M50(~has_M50_gaussian,:) = [];
GAUSSIAN_features.M100(~has_M100_gaussian,:) = [];
M50_gaussian(:,:,~has_M50_gaussian) = [];
M100_gaussian(:,:,~has_M100_gaussian) = [];

%%SAVE ALL FIGS

%Save data from STRF figures
save(strcat(figureFileStem, 'strf_features.mat'), 'STRF_features', 'GAUSSIAN_features');
save(strcat(figureFileStem, 'strf_data.mat'), 'STRF', 'M50_gaussian', 'M100_gaussian', '-v7.3');

%Save sources
save(M50_filename_peak, '-struct', 'sourceM50_peak');
save(M50_filename_zscore, '-struct', 'sourceM50_zscore');
save(M50_filename_lat, '-struct', 'sourceM50_lat');
%Save sources
save(M100_filename_peak, '-struct', 'sourceM100_peak');
save(M100_filename_zscore, '-struct', 'sourceM100_zscore');
save(M100_filename_lat, '-struct', 'sourceM100_lat');
%Save Gaussian sources
save(M50_filename_peak_GAUSSIAN, '-struct', 'sourceM50_peak_GAUSSIAN');
save(M50_filename_lat_GAUSSIAN, '-struct', 'sourceM50_lat_GAUSSIAN');
save(M50_filename_band_GAUSSIAN, '-struct', 'sourceM50_band_GAUSSIAN');
save(M100_filename_peak_GAUSSIAN, '-struct', 'sourceM100_peak_GAUSSIAN');
save(M100_filename_lat_GAUSSIAN, '-struct', 'sourceM100_lat_GAUSSIAN');
save(M100_filename_band_GAUSSIAN, '-struct', 'sourceM100_band_GAUSSIAN');
save(M50_filename_tempmod_GAUSSIAN, '-struct', 'sourceM50_tempMod_GAUSSIAN');
save(M100_filename_tempmod_GAUSSIAN, '-struct', 'sourceM100_tempMod_GAUSSIAN');


end


%% ===== Unique Filename Generator =====
function [filename, file_n] = filename_gen(fileStem, start_n)

if isempty(start_n)
    file_n = 1;
else
    file_n = start_n;
end

while exist(strcat(fileStem,num2str(file_n),'.mat'), 'file') == 2
     file_n = file_n + 1;
end
file_n = file_n + 1;
filename = strcat(fileStem,num2str(file_n),'.mat');
end

%% ===== STRF Generator =====
function strfMatrix = strfGen(compens1, iFreqEvents, OPTIONS, threshVals, events)
    
    %Create matrix structure
    strfMatrix = zeros(length(iFreqEvents), round(OPTIONS.strfLength/OPTIONS.binSize));
    
    for iEvent = 1:length(threshVals.TimeValues)
        
        timeZero = threshVals.TimeValues(iEvent);
        TimeWindow = [timeZero - (OPTIONS.strfLength/1000), timeZero];
        
        if TimeWindow(1,1) < events(strcmp({events.label},OPTIONS.allEventsName)).times(1) %Skip first event if its not at least 500ms from the start because not enough time to make an STRF
            continue;
        end
        
        for iFreq = 1:length(iFreqEvents)
            
            compens2 = threshVals.Values(iEvent); 
            finalCompens = compens1(iFreq)*compens2; % The value added to the STRF will be = to compens 2 (the z-score value of the associated activation event) * compens1
            
            allFreqTimes = events(iFreqEvents(iFreq)).times(events(iFreqEvents(iFreq)).times >= TimeWindow(1,1) & events(iFreqEvents(iFreq)).times < TimeWindow(1,2));
            allFreqTimes = allFreqTimes - timeZero + (OPTIONS.strfLength/1000); % Change the timepoint reference to be the STRF event
            allFreqTimes = ceil((allFreqTimes * 1000)/OPTIONS.binSize); % Change time resolution of STRF to binSize
            
            if~isempty(allFreqTimes)
                for ii = 1:length(allFreqTimes)
                    if allFreqTimes(ii) == 0 % In case the value was exactly on 0 and ceil didn't round up.
                        allFreqTimes(ii) = 1;
                    end
                    initialValue = strfMatrix(iFreq, allFreqTimes(ii));
                    strfMatrix(iFreq, allFreqTimes(ii)) = initialValue + finalCompens; % Apply compensation to the z-score value being added to STRF
                end
            end
        end
    end
    
    strfMatrix = strfMatrix ./ length(threshVals.TimeValues); % Convert STRF into average values
    
    % Post-processing:
    
    % Subtract mean of STRF baseline to centre baseline around zero
    baselineCorrect = strfMatrix(:,1:round(OPTIONS.strfBaseline/OPTIONS.binSize));
    meanMatrixCorrect = mean(baselineCorrect(:));
    strfMatrix = strfMatrix - meanMatrixCorrect;
    
    %Perform smoothing
    if OPTIONS.isMedfilt
        strfMatrix = medfilt2(strfMatrix, [3 3], 'zeros');
    end
    if OPTIONS.isSmooth
        strfMatrix = smoothdata(strfMatrix,1,'gaussian', 4);
        strfMatrix = smoothdata(strfMatrix,2,'gaussian', 4);
    end
end


%% ===== Gaussian-fit STRF Generator =====
function [M50_g, std_gaussianM50, M100_g, std_gaussianM100] = gaussianGen(strfMatrix, totalTTL, OPTIONS)

    %M50
    M50_g = zeros(size(strfMatrix,1),size(strfMatrix,2));
    [Ylim, Xlim] = size(strfMatrix);
    Yincrement = length(totalTTL)/Ylim; % length(totalTTL) being the number of frequencies there are
    yy = Yincrement:Yincrement:length(totalTTL);
    Xincrement = OPTIONS.strfLength/Xlim;
    xx = -OPTIONS.strfLength:Xincrement:-Xincrement;
    
    zz2 = strfMatrix(:,find(xx > -OPTIONS.gM50Window(2), 1):find(xx > -OPTIONS.gM50Window(1), 1));
    xx2 = 1:1:size(zz2,2);
    
    [fitresult, zfit] = fmgaussfit(xx2,yy,zz2);
    
    %Recreate a full 500ms matrix by assigning a zero value (in the
    %periphery of the gaussian) to the rest of the space.
    zeroVal = min(zfit(:));
    M50_g(:,1:find(xx > -OPTIONS.gM50Window(2), 1)-1)= deal(zeroVal);
    M50_g(:,find(xx > -OPTIONS.gM50Window(2), 1):find(xx > -OPTIONS.gM50Window(1), 1)) = zfit;
    M50_g(:,find(xx > -OPTIONS.gM50Window(1), 1)+1:end) = deal(zeroVal);
    std_gaussianM50 = fitresult(3:4)/sqrt(2);%standard deviation for X axis is in fitresult(3), Y axis in 4...needs to be corrected by sqrt 2 to give true standard dev
   
    %M100
    M100_g = zeros(size(strfMatrix,1),size(strfMatrix,2));
    zz2 = strfMatrix(:,find(xx > -OPTIONS.gM100Window(2), 1):find(xx > -OPTIONS.gM100Window(1), 1));
    xx2 = 1:1:size(zz2,2);
    
    [fitresult, zfit] = fmgaussfit(xx2,yy,zz2);
    
    %Recreate a full 500ms matrix by assigning a baseline value (in the
    %periphery of the gaussian) to the rest of the space.
    zeroVal = min(zfit(:));
    M100_g(:,1:find(xx > -OPTIONS.gM100Window(2), 1)-1) = deal(zeroVal);
    M100_g(:,find(xx > -OPTIONS.gM100Window(2), 1):find(xx > -OPTIONS.gM100Window(1), 1)) = zfit;
    M100_g(:,find(xx > -OPTIONS.gM100Window(1), 1)+1:end) = deal(zeroVal);
    std_gaussianM100 = fitresult(3:4)/sqrt(2); %standard deviation for X axis is in fitresult(3), Y axis in 4...needs to be corrected by sqrt 2 to give true standard dev
end


%% ===== STRF Feature Extractor =====
function [has_response, peakloc_col, peakloc_row, zscoreVal] = STRF_feat(subINT, meanMatrix, stdMatrix, strfLength, window, Xincrement)
    [M,I] = max(subINT(:));
    if M > 0
        has_response = true;
        [strfPeak(1), strfPeak(2)] = ind2sub(size(subINT), I);
        zscoreVal = (M - meanMatrix)/stdMatrix;
        peakloc_col = strfPeak(2) + (strfLength-window)/Xincrement;
        peakloc_row = strfPeak(1);
    else
        has_response = false;
        peakloc_col = [];
        peakloc_row = [];
        zscoreVal = [];
    end
end

%% ===== GAUSSIAN Feature Extractor =====
function [has_response, peakloc_col, peakloc_row, strfBand, strfTempMod] = GAUSSIAN_feat(subINT, std_gaussian, X, Y, OPTIONS, freqVals)
    
    [M,I] = max(subINT(:));
    
    if ~isempty(M) && M > 0
        has_response = true;
        
        % Peak
        [strfPeak(1), strfPeak(2)] = ind2sub(size(subINT), I);
        peakloc_col = strfPeak(2);
        peakloc_row = strfPeak(1);
        
        % Minimum in spectral domain
        strfMin = Y(peakloc_row) - (OPTIONS.baseDist*std_gaussian(2));
        if (strfMin < 1)
            strfMin = 1;
        end
            
        % Maximum in spectral domain
        strfMax = Y(peakloc_row) + (OPTIONS.baseDist*std_gaussian(2));
        if (strfMax > size(subINT,1))
            strfMax = size(subINT,1);
        end
            
        % Minimum in temporal domain    
        strfTMin = X(peakloc_col) - (OPTIONS.baseDist*std_gaussian(1));
            
        % Maximum in temporal domain
        strfTMax = X(peakloc_col) + (OPTIONS.baseDist*std_gaussian(1));
        
        % Bandwidth
        strfBand = log2(str2double(freqVals{ceil(strfMax)})) - log2(str2double(freqVals{ceil(strfMin)}));
        
        % Temporal modulation
        strfTempMod = 1/(2*(ceil(strfTMax) - ceil(strfTMin))/1000);
        
    else
        has_response = false;
    end

end

%% ===== GAUSSIAN FILTER FUNCTIONS =====
function [fitresult, zfit] = fmgaussfit(xx,yy,zz)
% FMGAUSSFIT Create/alter optimization OPTIONS structure.
%   [fitresult,..., rr] = fmgaussfit(xx,yy,zz) uses ZZ for the surface
%   height. XX and YY are vectors or matrices defining the x and y
%   components of a surface. If XX and YY are vectors, length(XX) = n and
%   length(YY) = m, where [m,n] = size(Z). In this case, the vertices of the
%   surface faces are (XX(j), YY(i), ZZ(i,j)) triples. To create XX and YY
%   matrices for arbitrary domains, use the meshgrid function. FMGAUSSFIT
%   uses the lsqcurvefit tool, and the OPTIMZATION TOOLBOX. The initial
%   guess for the gaussian is places at the maxima in the ZZ plane. The fit
%   is restricted to be in the span of XX and YY.
%   See:
%       http://en.wikipedia.org/wiki/Gaussian_function
%
%   Examples:
%     To fit a 2D gaussian:
%       [fitresult, zfit, fiterr, zerr, resnorm, rr] =
%       fmgaussfit(xx,yy,zz);
%   See also SURF, OMPTMSET, LSQCURVEFIT, NLPARCI, NLPREDCI.

%   Copyright 2013, Nathan Orloff. (Note, this function was shortened from
%   the original)

%% Condition the data
[xData, yData, zData] = prepareSurfaceData( xx, yy, zz );
xyData = {xData,yData};

%% Set up the startpoint
[amp, ind] = max(zData); % amp is the amplitude.
xo = xData(ind); % guess that it is at the maximum
yo = yData(ind); % guess that it is at the maximum
ang = 45; % angle in degrees.
sy = 1;
sx = 1;
zo = median(zData(:))-std(zData(:));
xmax = max(xData)+2;
ymax = max(yData)+2;
xmin = min(xData)-2;
ymin = min(yData)-2;

%% Set up fittype and options.
Lower = [0, 0, 0, 0, xmin, ymin, 0];
Upper = [Inf, 180, Inf, Inf, xmax, ymax, Inf]; % angles greater than 90 are redundant
StartPoint = [amp, ang, sx, sy, xo, yo, zo];%[amp, sx, sxy, sy, xo, yo, zo];

tols = 1e-16;
options = optimset('Algorithm','levenberg-marquardt',...
    'Display','off',...
    'MaxFunEvals',5e2,...
    'MaxIter',5e2,...
    'TolX',tols,...
    'TolFun',tols,...
    'TolCon',tols ,...
    'UseParallel','always');

%% perform the fitting
[fitresult,~,residual] = ...
    lsqcurvefit(@gaussian2D,StartPoint,xyData,zData,Lower,Upper,options);
zfit = gaussian2Duncert(fitresult,residual,xyData);
zfit = reshape(zfit,size(zz));

end

function zf = gaussian2Duncert(par,resid,xy)
% get the confidence intervals
J = guassian2DJacobian(par,xy);
[zf,~] = nlpredci(@gaussian2D,xy,par,resid,'Jacobian',J);
end

function z = gaussian2D(par,xy)
% compute 2D gaussian
z = par(7) + ...
    par(1)*exp(-(((xy{1}-par(5)).*cosd(par(2))+(xy{2}-par(6)).*sind(par(2)))./par(3)).^2-...
    ((-(xy{1}-par(5)).*sind(par(2))+(xy{2}-par(6)).*cosd(par(2)))./par(4)).^2);
end

function J = guassian2DJacobian(par,xy)
% compute the jacobian
x = xy{1}; y = xy{2};
J(:,1) = exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2);
J(:,2) = -par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*((2.*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))))./par(3).^2 - (2.*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))))./par(4).^2);
J(:,3) = (2.*par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2)./par(3)^3;
J(:,4) = (2.*par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2)./par(4)^3;
J(:,5) = par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*((2.*cosd(par(2)).*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))))./par(3).^2 - (2.*sind(par(2)).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))))./par(4).^2);
J(:,6) = par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*((2.*cosd(par(2)).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))))./par(4).^2 + (2.*sind(par(2)).*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))))./par(3).^2);
J(:,7) = ones(size(x));
end
