%% OPTICAT: Create Optimized ICA Training data 
% This function creates new "optimized" EEG dataset that can be used to
% train the ICA algorithmus to optimally model ocular ICA components during
% free viewing EEG experiments.  In the newly created, the sampling points 
% around certain experimental events (e.g. saccade onsets) are overweighted. 
% To do this, these samples will be repeatedly re-appended to the end of the 
% dataset. 
% This function is described and empircally evaluted in: 
%
% INPUTS
%  EEG    - EEG dataset (in EEGLAB format)
%           Can be epoched or continuous but must contain the event of
%           interest (e.g. 'saccade' events detected using EYE-EEG toolbox)
%
%  event2overweight - name of event in EEG.event. Data around this event
%           will be overweighted in the newly created dataset
%
%  timelim - [vector with two numbers] - time limits (in seconds) for the 
%            data epochs around the event2overweight to include. First value 
%            is before the event, second value is the time after the event.
%            Example: [-0.02 0.01] takes samples from -20 ms to +10 ms 
%            relative to the 'event2overweight'.
%
%  ow_proportion:
%     - this factor determines how many event samples are appended to the
%       end of the dataset, e.g. 0.1 = add  10% of length of orig. dataset
%                                1.0 = add 100% of length of orig. dataset
%                                (i.e. for a value of 1.0, the total number
%                                of samples in the dataset is doubled by 
%                                running this function)
% 
%  removemean - [boolean], remove mean value from each channel for the
%               appended, overweighted epochs?
%
% OUTPUTS
% - EEG with additional appended samples (=optimized ICA training dataset)
%
% Example call to this function:
% EEG_training = overweightEvents(EEG,'saccade',[-0.03 0.01],0.50,1)
%
% The function would create a dataset that consists of all the samples of
% the original EEG.data (in 2D) plus appended saccade samples. These
% appended samples comprise the the time interval from -30 ms to +10 ms
% relative to all 'saccade' events in the data. Samples around saccade
% events will be repeatedly re-appended to the end of the dataset until the
% dataset is 50% longer than before (ow_proportion). That is, the length of
% the newly created dataset is now 150% of its original length. The mean
% channel voltage will be removed from each appended short epoch
%
% Note 1: This function also subtracts the mean (baseline) of each channel 
% computed across the whole duration of each appended epoch. Removing the 
% epoch mean usually improves ICA decomposition.
%
% Note 2: The data should be optimally (high pass-) filtered for ICA 
% *before* running this function. See the reference paper:
%
% This procedure was proposed and empirically evaluated in: 
% Dimigen, O. (submitted). 
% Please cite this paper if you use this method.


function [EEG_out] = overweightevents(EEG,event2overweight,timelim,ow_proportion,removemean)

fprintf('\n\n%s(): Creating optimized ICA training dataset',mfilename);
fprintf('\n\n%s(): Creating dataset with event *%s* overweighted.\n', mfilename, event2overweight);
fprintf('%s(): Appended samples around event *%s* will make up %.f percent of original length of dataset.\n\n', mfilename, event2overweight, ow_proportion*100);
fprintf('%s(): This may take a moment...', mfilename);

disp('Inputs:')
event2overweight
[timelim]
ow_proportion
removemean


%% reshape original EEG.data to 2D format
EEG_tmp          = eeg_emptyset;
EEG_tmp.data     = EEG.data(:,:); % reshape data to 2D
EEG_tmp.pnts     = size(EEG_tmp.data,2);
EEG_tmp.srate    = EEG.srate;
EEG_tmp.chanlocs = EEG.chanlocs;
EEG_tmp.event    = EEG.event;
EEG_tmp.nbchan   = size(EEG_tmp.data,1);
EEG_tmp          = eeg_checkset(EEG_tmp);

%% how overweight event samples to add?
npoints    = size(EEG_tmp.data,2); % number of EEG samples
nsacpoints = round(ow_proportion * npoints); % desired number of (redundant) saccade samples

%% create event-locked epochs to overweight
sac = pop_epoch(EEG,{event2overweight},timelim);
if removemean
    sac = pop_rmbase(sac,[]); % baseline is subtracted across whole epoch
end
%sac = applytochannels(sac,1:NCHANS_EEG,'pop_rmbase( EEG,[]);');

%% overweight (=copy & re-append) overweight-event-locked epochs
tmpsacdata = sac.data(:,:); % 2D
tmpsacdata = repmat(tmpsacdata,1, ceil(nsacpoints/size(tmpsacdata,2)));
tmpsacdata = tmpsacdata(:,1:nsacpoints);
clear sac

%% create EEG dataset for overweight-event-locked epochs
EEG_sac          = eeg_emptyset;
EEG_sac.data     = tmpsacdata;
EEG_sac.pnts     = size(EEG_sac.data,2);
EEG_sac.srate    = EEG.srate;
EEG_sac.chanlocs = EEG.chanlocs;
EEG_sac.nbchan   = size(EEG_sac.data,1);
EEG_sac.event    = []; % overweighted events not in EEG.event structure
EEG_sac          = eeg_checkset(EEG_sac);

%% merge original EEG dastset with overweighted event-locked epochs
EEG_out = pop_mergeset(EEG_tmp,EEG_sac,1);

fprintf('%s(): Done.', mfilename);