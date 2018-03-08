% pop_importeyetracker() - import and synchronize simultaneously recorded
%           eye tracking data and store it as additional data channels in
%           the EEG structure. Loads the MATLAB copy of the raw ET data
%           generated by pop_parsesmi() or pop_parseeyelink().
%           Synchronizes ET and EEG based on common events present in both
%           recordings. The eye tracking data will be resampled/linearly
%           interpolated to match the sampling frequency of the EEG.
%           If available, eye movement events (blink/saccade/fixation
%           detected online) can be imported from the ET raw data (for
%           Eyelink only).
%
% Usage:
%   >>  EEG = pop_importeyetracker(EEG, matFileToLoad, startEndEvent, ...
%   importColumns, newLabels, importEyeEvents, doRegression, ...
%   filterEyetrack, plotFig)
%
% Inputs:
%   EEG	            - EEG dataset (continuous data)
%                     for which there is a parallel ET recording
%   matFileToLoad   - <string> file path to eyetracker mat-file
%                     generated by function parsesmi or parseeyelink
%   startEndEvent   - <vector with two values> event codes shared in the
%                     EEG and ET recording for synchronization. A vector
%                     with two numerical values is expected:
%                     [startevent, endevent], e.g. [102 202]
%   importColumns	- <double> array of indices of eyetracker columns to
%                     import, e.g., [3 5 7]
%   newLabels	    - <cell> of <string> of length(importColumns) which
%                     containts the channel labels of the newly
%                     imported ET data columns,
%                     e.g. {'gaze_x', 'gaze_y', 'pupilsize'}
%   importEyeEvents	- <boolean> [0|1] set TRUE if you want to import
%                     eye movement events (e.g. saccades) that were
%                     detected online by the ET recording software
%                     (works for Eyelink only). Imported events will be
%                     stored in the EEG.events structure.
%   doRegression    - <boolean> [0|1]
%                     FALSE: interpolate data only based on the latency of
%                     the start-event and end-event
%                     TRUE: estimate latency of start-event and end-event
%                     based on linear regression using all events that are
%                     shared by the EEG and ET recording. For GUI use, the
%                     default setting for doRegression is TRUE.
%   filterEyetrack  - <boolean> [0|1] {default: 0}
%                     FALSE: resample the eye tracking data to the sampling
%                     frequency of the EEG (using linear interpolation),
%                     but do not apply any filters to the eye tracking data
%                     TRUE: filter the eye tracking data at the Nyquist
%                     frequency to prevent aliasing/image artifacts that
%                     may result from downsampling or upsampling the ET
%                     data to the sampling frequency of the EEG. This
%                     function is NOT YET IMPLEMENTED so the setting has no
%                     effect at the moment.
%
%   plotFig	        - <boolean> [0|1] set TRUE to plot figure showing the
%                     quality of synchronisation (recommended)
%
%   searchRadius    - radius (in samples) around EEG events to search for
%                     matching (shared) events of the same type in the ET.
%                     This value is important to identify instances of 
%                     shared events that were transmitted between the 
%                     startEvent and endEvent.
%                     The algorithm will search within plusminus
%                     'searchRadius' samples around each EEG event for an
%                     ET event of the same type (e.g. '123').
%                     Unless you have trouble with synchronization, or use
%                     a very low (e.g. 100 Hz) or very high (e.g. 1000 Hz)
%                     sampling rate for the EEG, we recommend sticking to 
%                     the default value of 4.
%                     However, it may be a good idea to increase the
%                     searchRadius in case of very high sampling rates
%                     (e.g. 2000 Hz) or in the case of errors that can happen
%                     if multiple events of the same type (e.g. '123') were
%                     send in immediate succession (that is, less than
%                     2*searchRadius+1 samples apart) during the experiment.
%                     [DEFAULT: 4 samples]
%
%
% Outputs:
%   EEG             - EEG dataset (with synchronized eye tracking data)
%
% See also:
%   pop_parsesmi, pop_parseeyelink, pop_parsetobii, synchronize
%
% Example: An call of pop_importeyetracker() might look like this:
% >> EEG = pop_importeyetracker(EEG,'C:\\eyedata.mat',[103 203],[3 7],...
%          {'Gaze_x' 'Pupil_Dia_(mm)'},0,1,0,5,1)
%
% Explanation: load the parsed ET data stored in "eyedata.mat" that was
% generated by pop_parsesmi, pop_parseeyelink, or pop_parsetobii. Use the 
% first event of type "103" and the last event of type "203" as the startEvent
% and endEvent for synchronization, respectively. Import the third and
% seventh numerical data column from the ET raw data. These columns contain
% the horizontal gaze position and pupil diameter and are named accordingly.
% Do not import eye movement events from the raw data. Do a linear
% regression on all shared events to optimize the latency estimation for
% the startevent and endevent. Do not filter the eye tracking data.
% Plot a figure that informs about the precision of synchronization. Set
% the searchRadius to identify matching triggers to 5 samples.
%
% Author: ur & od
% Copyright (C) 2009-2017 Olaf Dimigen & Ulrich Reinacher, HU Berlin
% olaf.dimigen@hu-berlin.de

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, 51 Franklin Street, Boston, MA 02110-1301, USA


function [EEG, com] = pop_importeyetracker(EEG,matFileToLoad,startEndEvent,importColumns,channelLabels,importEyeEvents,doRegression,filterEyetrack,plotFig,searchRadius)

%% display help if no arguments
if nargin < 1
    help(mfilename);
    return;
end

%% initialize history command string
com = '';

%% test whether a continuous EEG dataset is already loaded
if isempty(EEG.data)
    warndlg2('The corresponding EEG file must be loaded first to synchronize it with eye tracking data!',[mfilename '(): Dataset is empty!'])
    return
elseif size(EEG.data,3) > 1 % data is epoched
    warndlg2('For synchronization with eye tracking data, the EEG must be continuous (not yet epoched)!',[mfilename '(): Dataset is epoched!'])
    return
end

%% look for boundary events in EEG data
ixBnd = strcmp('boundary',{EEG.event.type});
if ixBnd
    if isfield(EEG.event,'duration') && any( [EEG.event( ixBnd ).duration] > 1)
        warndlg2(sprintf('You have event(s) of type ''boundary'' that represent removed samples (duration > 1).\n This will cause problems if those boundaries are in the synchronization range.'),[mfilename '(): Boundaries in Continuous Data'])
    else
        warndlg2(sprintf('You have event(s) of type ''boundary'' without duration information. Those might represent removed samples (duration > 1)\n and will cause problems if they in the synchronization range.'),[mfilename '(): Boundaries in Continuous Data'])
    end
end

%% display second GUI dialog (set synch parameters)?
popDialog2 = false;

try
    %% load parsed ET data
    if ~exist('matFileToLoad','var')
        % if matFileToLoad is not provided, pop up dialogue
        matFileToLoad = dlg_loadeyetracker(mfilename);
    end;
    
    %% load parsed ET structure with data and information
    if exist(matFileToLoad,'file')
        fprintf('\n%s(): Loading %s...',mfilename,matFileToLoad)
        ET = load(matFileToLoad);
    else
        error('\n%s: Did not find a file named %s. Please check filepath or filename.',mfilename,matFileToLoad);
    end
    
    %% are ET events (triggers or messages containing keyword) present in ET.event?
    if ~isfield(ET,'event')
        ET.event = [];
    end
    if isempty(ET.event)
        fprintf('\n%s(): Found no trigger/event information in your parsed eye tracking data (ET.event)!\n',mfilename)
        fprintf('\nSee the help of the pop_parse >your Eyetracker< function and read the plugin description.\n')
        return
    end
    
    %% are eye movement events in the data (detected online by ET system)? [Eyelink only]
    hasEyeEvents = isfield(ET,'eyeevent');
    if ~exist('importEyeEvents','var')
        importEyeEvents = 0;
    end
    if ~hasEyeEvents && importEyeEvents
        fprintf('%s(): Found no eye movement events (blinks/saccades/fixations) in the raw data!',mfilename)
    end
    
    %% retrieve EEG trigger from EEG.event
    % check EEG.event.type for triggers present in the EEG recording
    % sometimes EEG triggers are represented as a string ("S 123") rather than
    % an integer value [123], e.g. for data recorded with BrainProducts amplifiers
    
    if any(cellfun(@isstr,{EEG.event(:).type})) % strings in EEG.event.type
        % find consecutive digits in event codes, e.g. 'S 123', 'R234',
        tmp = regexp({EEG.event.type},'\d+','match')';
        % find events without any numerical value, e.g. 'Start', 'Boundary'
        noTrigValues = find(cellfun(@isempty,tmp(:)));
        eegEvents(:,2) = [EEG.event.latency]';
        
        % clear empty strings and their latencies
        tmp(noTrigValues) = [];
        eegEvents(noTrigValues,:) = [];
        eegEvents(:,1) = str2double([tmp{:}]');
        
        clear tmp noTrigValues
    else % no strings in EEG.event.type
        eegEvents(:,1) = [EEG.event.type]';
        eegEvents(:,2) = [EEG.event.latency]';
    end
    % eegEvents: [type, latency, latencyInSyncRange]
    
    %% find events shared by EEG and eye tracker
    eegTypes    = unique(eegEvents(:,1));
    etTypes     = unique(ET.event(:,2));
    sharedTypes = intersect(eegTypes,etTypes);
    if isempty(sharedTypes)
        fprintf('\n%s(): There are no shared events that occur both in the EEG and the eye track!\n',mfilename)
        return
    end
    
    %% prepare and pop dialogue
    % - to select start-event and end-event for synchronization
    % - to specify which ET columns to import
    % - to define names (channel labels) of newly imported ET channels
    if exist('startEndEvent','var')
        [syncTypesPresent,selectInPullDown] = ismember(startEndEvent,sharedTypes);
        
        if any(~syncTypesPresent)
            error('\n%s(): Did not find events of the specified type [%s] in both ET and EEG data!\nPlease check your raw data.',mfilename,num2str(startEndEvent(~syncTypesPresent)));
        end
    else
        selectInPullDown = [1 1]; % default behaviour in pulldown menu
        popDialog2 = true;
    end
    
    nrColumns = size(ET.data,2);
    
    if ~exist('importColumns','var')
        importColumns = 1:nrColumns;
        popDialog2 = true;
    end
    
    if ~exist('channelLabels','var')
        channelLabels = ET.colheader(importColumns);
    end;

    if ~exist('doRegression','var')
        doRegression = true;
    end    
    
    if ~exist('plotFig','var')
        popDialog2 = true;
    end

    if ~exist('searchRadius','var')
        searchRadius = 4; % default search radius is 4 samples around EEG trigger
    end    
    
    if popDialog2
        % display the number a specific event type in the ET and EEG data
        nrEeg = hist(eegEvents(:,1),eegTypes);
        nrEeg = nrEeg(find(nrEeg)); % protect against hist() behaviour for 1 category
        nrEeg = nrEeg(ismember(eegTypes, sharedTypes));
        nrEt = hist(ET.event(:,2),etTypes);
        nrEt = nrEt(find(nrEt)); % protect against hist() behaviour for 1 category
        nrEt = nrEt(ismember(etTypes, sharedTypes));
        
        % for GUI dialogue, display example data from each ET data column
        % For this, take first non-zero element of each data column.
        % Otherwise zero.
        [dummy,dummy,exampleData] = arrayfun(@(x) find(ET.data(:,x),1),1:nrColumns,'UniformOutput',false);
        exampleData(cellfun(@(x) isempty(x),exampleData)) = {0};
        
        % pop window to get user input
        [startEndEvent,importColumns,channelLabels,importEyeEvents,plotFig] = ...
            dlg_syncsettings(mfilename, sharedTypes, nrEeg, nrEt, selectInPullDown, nrColumns, importColumns, channelLabels, exampleData, hasEyeEvents, importEyeEvents);
        clear nrEeg nr* Et dummy
    end
    
    % default for GUI use: set filterEyetrack to false (still experimental)
    if ~exist('filterEyetrack','var')
        filterEyetrack = false;
    end
    
    %% synchronize ET and EEG based on start-event and end-event
    fprintf('\n%s(): Synchronizing EEG and eye tracking data...',mfilename);
    [ET,eegEvents,syncErrorTable] = synchronize(ET, startEndEvent, eegEvents, EEG.srate, EEG.pnts, doRegression, filterEyetrack, plotFig, searchRadius);
    
    %% store info about sync quality with dataset (EEG.etc)
    EEG.etc.eyetracker_syncquality = syncErrorTable;
    
    %% add ET data and new channel labels to EEG structure
    % test for illegal characters in channel names incompatible with EEGLAB
    changedName = false;
    for n = 1:length(channelLabels)
        if strfind(channelLabels{n},' ') | strfind(channelLabels{n},'[') | strfind(channelLabels{n},']') % spaces, brackets
            changedName = true;
        end
        channelLabels(n) = strrep(channelLabels(n),' ','-'); % causes various problems % BUGFIX 09-2015: use to be replaced with '_'
        channelLabels(n) = strrep(channelLabels(n),'_','-'); % causes problem with eeg_decodechans() % BUGFIX 09-2015: use to be replaced with '_'
        channelLabels(n) = strrep(channelLabels(n),']',')'); % causes problem with eeg_decodechans()
        channelLabels(n) = strrep(channelLabels(n),'[','('); % causes problem with eeg_decodechans()
        
    end
    if changedName
        fprintf('\n%s(): Replaced characters in ET channel names for compatibility with EEGLAB\n',mfilename)
    end
    for col=1:length(importColumns)
        EEG.chanlocs(EEG.nbchan +col).labels = channelLabels{col};
        EEG.chanlocs(EEG.nbchan +col).ref  = '';
        EEG.chanlocs(EEG.nbchan +col).type = 'EYE';
    end
    
    % add new ET channels
    if ~isempty(ET.syncdata)
        sampleFirstEvent = eegEvents(1,2);
        sampleLastEvent  = eegEvents(end,2);
        
        % set all ET data outside the synchronized range to zero
        EEG.data = [EEG.data; zeros(length(importColumns),EEG.pnts)];
        EEG.data(EEG.nbchan+1:end, sampleFirstEvent:sampleLastEvent) = ET.syncdata(:,importColumns)';
        EEG.nbchan = size(EEG.data,1);
    end
    
    %% import eye movement events from raw data (Eyelink only)
    if hasEyeEvents && importEyeEvents
        fprintf('\n%s(): Importing eye movement events from eyetracker...',mfilename);
        eventTypes = fieldnames(ET.eyeevent)';
        importedEyeTypes = {};
        n_imported       = [];
        % loop tru event types (saccades,fixations,blinks)
        for evtt = 1:length(eventTypes)
            eventType = eventTypes{evtt};
            eyeCodes  = unique(ET.eyeevent.(eventType).eye);
            % loop tru eyes
            for thisEye = 1:length(eyeCodes)
                % get index for left/right eye events
                ix_thisEye = strcmp(eyeCodes(thisEye),cellstr(ET.eyeevent.(eventType).eye));
                eyefieldnames = ET.eyeevent.(eventType).colheader;
                eyeType = [eyeCodes(thisEye) '_' eventType(1:end-1)];
                % add eye events to EEG.event
                EEG = addevents(EEG,ET.eyeevent.(eventType).data(ix_thisEye,:),eyefieldnames,eyeType);
                % remember which/how many events were imported
                importedEyeTypes = [importedEyeTypes; cellstr(eyeType)];
                n_imported       = [n_imported; length(ix_thisEye)];
            end
        end
        fprintf('\nThe following eye movement event types were imported:\n%s\n',vararg2str(importedEyeTypes))
        for j = 1:length(importedEyeTypes)
            fprintf('\n %i events of type %s:',n_imported(j),importedEyeTypes{j});
        end
        fprintf('\n');
        
        clear importedEyeTypes
    end
    
    EEG = eeg_checkset(EEG,'eventconsistency');
    
    %% store infos about 'other' ET message lines in the EEG structure
    % other messages are ET messages that are not 'eyeevents'
    % (see above) and not special keyword-messages (used for
    % synchronization), but which may contain some aspects of the 
    % experimental design
    if isfield(ET,'othermessages')
        EEG.etc.eyetracker_othermessages = ET.othermessages;
        fprintf('\nImporting infos about other eyetracking messages into EEG.etc\n');
    end
    
catch err
    if (strcmp(err.identifier,'MATLAB:unassignedOutputs'))
        % disp('User decided to cancel the selection of input.')
        return
    else
        rethrow(err);
    end
end

%% return string for command history
allArgs = vararg2str({matFileToLoad, startEndEvent, importColumns, channelLabels, importEyeEvents, doRegression, filterEyetrack, plotFig, searchRadius});
com = sprintf('EEG = %s(EEG,%s)',mfilename,allArgs);
fprintf('\nSynchronization completed.\n')