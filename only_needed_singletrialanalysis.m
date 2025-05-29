
%%%%%%%%%%%%%%%%%%%%%%% DATA Loading %%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/storage/share/matlab/labbox/'))

% Define FileBase and FileIn correctly
FileBase = '/storage3/manus/Openephys Recordings/';
RatID = 'MP02';
TrialID = [RatID '-20251403-trial19'];
FileTrial = [FileBase TrialID];
FileIn = [FileTrial '/Record Node 102/experiment2/recording1/continuous/Acquisition_Board-100.Rhythm Data/continuous.dat'];

% Load LFP data
lfp = LoadBinary(FileIn, 1:8, 8);

% Load parameters
xmlFile = [FileBase 'Stats.xml'];
par = LoadXml(xmlFile);
lfpSamplingRate = par.SampleRate; % Sampling rate in Hz

% Extract rotary encoder signals (Channels 2 & 3)
encoder_signal_2 = lfp(2, :); 
encoder_signal_3 = lfp(3, :); 


% Average both signals for better accuracy
encoder_signal = (encoder_signal_2 + encoder_signal_3) / 2;

%%%%%%%%%%%%%%%%%%%%%%% Speed Computation %%%%%%%%%%%%%%%%%%%%%%

% Improved threshold calculation (using 70% of signal range)
min_signal = min(encoder_signal);
max_signal = max(encoder_signal);
threshold = min_signal + 0.7 * (max_signal - min_signal);

% Detect rising edges (pulses)                                  <==== This approach  picks the pulse onsets alternating either on one channel or on the
% another. Not sure it is correct!!
pulse_indices = find(diff(encoder_signal > threshold) == 1); % Get pulse times

% Compute time differences (time between pulses)
pulse_times = pulse_indices / lfpSamplingRate; % Convert indices to time in seconds
pulse_intervals = diff(pulse_times); % Time difference between consecutive pulses

% Compute speed (in cm/s)
wheel_radius_cm = 18; % Corrected wheel radius
circumference_cm = 2 * pi * wheel_radius_cm;

%Corrected by Evgeny
rotations_per_sec = (1/600) ./ pulse_intervals;
speed_cm_per_sec = (rotations_per_sec * circumference_cm);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%====== Part added by Evgeny ===========================================================================%
% %Evgeny's way to compute speed (cm/sec) as a DOUBLE-CHECK
% %For now i take the detected pulse onset times from Manus to check the speed calculation
% %Compute inter-pulse time intervals (sec)
% InterpulseIntervals_sec = diff(pulse_times);
% %Number of pulses per revolution (from the rotary encoder spec)
% nPulsesPerRevolution = 600;
% %Angle increment per pulse (degrees)
% AngularIncrement = 360/nPulsesPerRevolution;
% %Compute the instanteneous angular speed (degree/sec)
% AngSpeed = AngularIncrement ./ InterpulseIntervals_sec;
% %Radius of the wheel (cm) - SET BY MANUS ABOVE TO 18 cm!
% %Compute the instanteneous linear speed (sm/sec)
% Speed = AngSpeed * (2*pi*wheel_radius_cm/360);
% %Smooth the speed by moving average kernel
% Speed_smoothed = smooth(Speed,20);


%Create a time vector for the TTL signal (sec)
time_encoder = [0:length(encoder_signal_2)-1]/lfpSamplingRate;

%Create a time vector for speed (sec)
time_speed = pulse_times(1:end-1);

%Detect and discard speed artifact due to a technical issue in the TTL pulse signal
%NOTE: these artifacts seem to be due to some technical issue of the signal from the rotary
%encoder, when the inter-pulse time interval becomes all of a sudden drops.

%Compute instanteneous acceleration (cm/sec^2) using the raw (non-smoothed) speed vector
acc = diff(speed_cm_per_sec) ./ pulse_intervals(1:end-1);
%Create a time vector for acceleration (sec)
time_acc = time_speed(1:end-1);

%Detect NEGATIVE acceleration peaks (of any length)
clear AccPeak_neg
%Threshold for the animal acceleration (by absolute value)
ThrdAcc_neg = -10000;
%Binary map of suprathreshold acceleration (cm/sec^2)
binmap = acc < ThrdAcc_neg;
%Calculate start/end indices of the periods of continuous suprathreshold acceleration 
out = contiguous(binmap,1);
if ~isempty(out) && numel(out) >= 2 && ~isempty(out{1,2})     per = out{1,2}; else     per = []; end 
per_sec = time_acc(per);
%In each period peak the sample with the largest amplitude
for n=1:size(per,1)
    ind = per(n,1):per(n,2);
    [~,min_ind] = min(acc(ind));
    AccPeak_neg(n) = ind(min_ind);
end
clear out ind min_ind n binmap per per_sec


%Detect POSITIVE acceleration peaks (of any length)
clear AccPeak_pos
%Threshold for the animal acceleration (by absolute value)
ThrdAcc_pos = 3000;
%Binary map of suprathreshold acceleration (cm/sec^2)
binmap = acc > ThrdAcc_pos;
%Calculate start/end indices of the periods of continuous suprathreshold acceleration 
out = contiguous(binmap,1);
if ~isempty(out) && numel(out) >= 2 && ~isempty(out{1,2})     per = out{1,2}; else     per = []; end 
per_sec = time_acc(per);
%In each period peak the sample with the largest amplitude
for n=1:size(per,1)
    ind = per(n,1):per(n,2);
    [~,min_ind] = max(acc(ind));
    AccPeak_pos(n) = ind(min_ind);
end
clear out ind min_ind n binmap per per_sec


%Time periods which include one or more artifact acceleration peaks within
%a short time interval (0.25 sec) are defined as bad periods.
%Merge positive and negative peaks together
AccPeak = union(AccPeak_neg, AccPeak_pos)';
% AccPeak_sec = time_acc(AccPeak)';
%Compute the interpeak intervals (sec)
InterPeakIntervals = diff(time_acc(AccPeak))';
%Binary map of all interpeak intervals which are <0.25 sec
binmap = InterPeakIntervals < 0.25;
%Find contiguous sequences of suprathreshold (short) interpeak intervals
out = contiguous(binmap,1);
if ~isempty(out) && numel(out) >= 2 && ~isempty(out{1,2})     per = out{1,2}; else     per = []; end 
%Periods (start/end) of the speed samples to be discarded
SpeedSamplePeriods2Discard = AccPeak([per(:,1) per(:,2)+1]); 
%Speed samples to be discarded
SpeedSamples2Discard = [];
for n=1:size(SpeedSamplePeriods2Discard,1)
    ind = SpeedSamplePeriods2Discard(n,1) : SpeedSamplePeriods2Discard(n,2);
    SpeedSamples2Discard = [SpeedSamples2Discard ind];    
end
%Add to SpeedSamples2Discard also single peaks which do not belong to any contiguous sequences 
ind = ~WithinRanges(AccPeak, SpeedSamplePeriods2Discard);
SingleAccPeaks2Discard = AccPeak(ind) + 1; %+1 to properly refer to speed samples
SpeedSamples2Discard = sort([SpeedSamples2Discard SingleAccPeaks2Discard']);
clear ind binmap out per n 


%Discard the speed samples in the bad periods (replace with NaN)
speed_cln = speed_cm_per_sec;
speed_cln(SpeedSamples2Discard) = NaN; 


%Interpolate the speed samples using a new time vector with a constant dt
%NOTE: interpolation is only done during the RUN periods!!!!!
SpeedSamplingRate = 50;
time_speed2 = 0 : (1/SpeedSamplingRate) : time_encoder(end);
%Interpolate the speed samples using a new time vector with a constant dt
%IMPORTANT: The very first available speed sample happens at time >0, 
%when the first TTL pulse from the rotary encoder happened.
%Hence, the interpolation must be done between the very first and last
%speed samples to avoid edge artifacts.
%TO DO: interpolation must only be done during the RUN periods, when speed was defined.
% Set interpolation rate (Hz)
SpeedSamplingRate = 50;

% Time vector for interpolation (uniform sampling)
time_speed2 = 0 : (1/SpeedSamplingRate) : time_encoder(end);

% Initialize interpolated speed with NaNs
speed_interp = NaN(size(time_speed2));

% Determine valid interpolation window (avoid extrapolation)
FirstSpeedSampleTime = time_speed(1);
LastSpeedSampleTime  = time_speed(end);

% Optional: trim padding to avoid unstable edges
padding_sec = 1;  % seconds to trim from both ends
valid_interp_mask = time_speed2 >= (FirstSpeedSampleTime + padding_sec) & ...
                    time_speed2 <= (LastSpeedSampleTime - padding_sec);

% Perform shape-preserving interpolation with extrapolation set to NaN
interpspeed_insert = interp1(time_speed, speed_cln, time_speed2(valid_interp_mask), 'pchip', NaN);

% Insert interpolated values
speed_interp(valid_interp_mask) = interpspeed_insert;
% speed_interp(isnan(speed_interp)) = min(speed_cln);

% Clip interpolated values to min/max of clean speed (safety clamp)
%speed_interp = max(min(speed_interp, nanmax(speed_cln)), nanmin(speed_cln));


% ==== Gaussian smoothing after interpolation ====
window_size = 40;         % still 40 samples wide
gauss_kernel = gausswin(window_size);
gauss_kernel = gauss_kernel / sum(gauss_kernel);



% Apply smoothing
speed_smth = conv(speed_interp, gauss_kernel, 'same');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Define inactivity periods
min_pause_duration = 5; % Define inactivity threshold (e.g., 2 seconds)
[inactivity_start, inactivity_end] = detect_rotary_encoder_activity(pulse_times, pulse_indices, lfpSamplingRate, encoder_signal, min_pause_duration);




%%%%%% LOW threshold start and end

% Define the new low-speed threshold
low_speed_threshold = 8; 

% Identify indices where speed is above threshold
above_threshold = speed_cln > low_speed_threshold;

% Find transitions (start and end of running periods)
start_indices = find(diff([0, above_threshold]) == 1); % Rising edge (start of running)
end_indices = find(diff([above_threshold, 0]) == -1); % Falling edge (end of running)

% Extract corresponding times
running_periods = [pulse_times(start_indices)', pulse_times(end_indices)'];

% Display the matrix
disp('Start and End times of real running behavior:');
disp(running_periods);


%%%%%%%%%%%% INACTIVE PERIODS %%%%%%%%%%%%%%%%%%%%%%%%

% Apply Gaussian smoothing only to active periods
active_mask = true(size(speed_cln)); % Assume all values are active initially

% Mark inactive indices
for i = 1:length(inactivity_start)
    inactive_indices = time_speed2 >= inactivity_start(i) & time_speed2 <= inactivity_end(i);
    active_mask(inactive_indices) = false; % Mark inactive periods
end

% Extract only active speed values
active_speeds = speed_cln(active_mask);  

% Apply Gaussian smoothing to active speeds only
smoothed_active_speed = conv(active_speeds, gauss_kernel, 'same');

% Ensure inactive periods remain exactly zero
speed_smth(~active_mask) = 0;

%%%%%%% Group Running Periods if Separation < 2 seconds %%%%%%%
min_gap = 2; % Minimum gap to merge running periods (in seconds)
merged_running_periods = [];
current_start = running_periods(1,1);
current_end = running_periods(1,2);

for i = 2:size(running_periods, 1)
    if running_periods(i,1) - current_end < min_gap
        % Merge this period with the current one
        current_end = running_periods(i,2);
    else
        % Save the merged period and start a new one
        merged_running_periods = [merged_running_periods; current_start, current_end];
        current_start = running_periods(i,1);
        current_end = running_periods(i,2);
    end
end
% Add the final period
merged_running_periods = [merged_running_periods; current_start, current_end];

%%%%%%% Remove running episodes shorter than 5 seconds %%%%%%%
min_running_duration = 5; % seconds
durations = merged_running_periods(:,2) - merged_running_periods(:,1);
keep_idx = durations >= min_running_duration;
merged_running_periods = merged_running_periods(keep_idx, :);

%%%%%%% Detect Frequent Direction Changes and Remove Periods %%%%%%%

% Define parameters
direction_change_threshold = 5; % More than 5 changes
window_duration = 5; % In seconds  

% Identify direction at each pulse
direction_mask = encoder_signal_2(pulse_indices(1:end-1)) > encoder_signal_3(pulse_indices(1:end-1));

 
direction_interp = NaN(size(valid_interp_mask));

% Interpolate with extrapolation to fill the beginning
direction_interp(valid_interp_mask) = interp1(time_speed, double(direction_mask), ...
                                              time_speed2(valid_interp_mask), ...
                                              'previous', NaN);
direction_interp(direction_interp == 0) = -1;

direction_interp = conv(direction_interp, gauss_kernel, 'same');
direction_interp(direction_interp >= 0) = 1;
direction_interp(direction_interp < 0) = -1;

speed_with_direction = direction_interp.*speed_smth;

% Detect when direction changes
direction_changes = find(diff(direction_mask) ~= 0);
direction_change_times = pulse_times(direction_changes); % Convert indices to time

% Find periods with too many direction changes
remove_periods = [];
for i = 1:length(direction_change_times)
    % Find changes within the window duration
    count_changes = sum(direction_change_times >= direction_change_times(i) & direction_change_times <= direction_change_times(i) + window_duration);
    
    if count_changes > direction_change_threshold
        % Mark this window for removal
        remove_periods = [remove_periods; direction_change_times(i), direction_change_times(i) + window_duration];
    end
end

% Merge overlapping remove periods
if ~isempty(remove_periods)
    merged_remove_periods = [remove_periods(1,1), remove_periods(1,2)];
    for i = 2:size(remove_periods, 1)
        if remove_periods(i,1) <= merged_remove_periods(end,2) % Overlapping
            merged_remove_periods(end,2) = max(merged_remove_periods(end,2), remove_periods(i,2));
        else
            merged_remove_periods = [merged_remove_periods; remove_periods(i,1), remove_periods(i,2)];
        end
    end
else
    merged_remove_periods = [];
end

%%%%%%% Trim Running Periods Instead of Removing Completely %%%%%%%
filtered_running_periods = [];
for i = 1:size(merged_running_periods, 1)
    current_start = merged_running_periods(i,1);
    current_end = merged_running_periods(i,2);
    
    % Check for overlapping filtered-out segments
    for j = 1:size(merged_remove_periods, 1)
        remove_start = merged_remove_periods(j,1);
        remove_end = merged_remove_periods(j,2);
        
        if remove_end <= current_start || remove_start >= current_end
            continue; % No overlap, skip
        elseif remove_start <= current_start && remove_end >= current_end
            current_start = NaN; % Whole period removed
            break;
        elseif remove_start > current_start && remove_end < current_end
            % Split into two valid segments
            filtered_running_periods = [filtered_running_periods; current_start, remove_start];
            current_start = remove_end;
        elseif remove_start <= current_start
            current_start = remove_end;
        elseif remove_end >= current_end
            current_end = remove_start;
        end
    end
    
    % Add the remaining valid portion if still valid
    if ~isnan(current_start) && current_start < current_end
        filtered_running_periods = [filtered_running_periods; current_start, current_end];
    end
end

%%%%%%% Final cleanup: Remove fragments shorter than 5 seconds after trimming %%%%%%%
min_running_duration = 5; % seconds
if ~isempty(filtered_running_periods)
    durations = filtered_running_periods(:,2) - filtered_running_periods(:,1);
    keep_idx = durations >= min_running_duration;
    filtered_running_periods = filtered_running_periods(keep_idx, :);
end


% Remove data from the first second
valid_idx = time_speed2 > 1;

%%%%% Rewards %%%%%%%%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Extract channel 8 (reward delivery TTL pulses)
reward_signal = lfp(8, :);

% Define a threshold for detection (assuming binary TTL pulses)
reward_threshold = max(reward_signal) * 0.5; % Adjust if necessary

% Detect rising edges (TTL pulses indicating rewards)
reward_indices = find(diff(reward_signal > reward_threshold) == 1);

% Convert indices to time (in seconds)
reward_timestamps = reward_indices / lfpSamplingRate;

% add the reward times -5s to the indices
for rew_itr = 1:length(reward_timestamps)
   
    t = reward_timestamps(rew_itr);
    if ~any(t >= filtered_running_periods(:,1) & t <= filtered_running_periods(:,2))
        filtered_running_periods = [filtered_running_periods; t-5, t];
    end
    
end

% sort
filtered_running_periods = sortrows(filtered_running_periods,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Excel %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% COMPUTE EPISODE METRICS %%%%%%%%%%%%%%%%%%%%%%%

episode_starts = filtered_running_periods(:,1);
episode_ends = filtered_running_periods(:,2);
episode_durations = episode_ends - episode_starts;

% Preallocate episode metrics
episode_paths = zeros(size(episode_durations));
episode_directions = zeros(size(episode_durations));
episode_mean_speeds = zeros(size(episode_durations));

for i = 1:length(episode_starts)
    mask = time_speed2 >= episode_starts(i) & time_speed2 <= episode_ends(i);
    episode_speeds = speed_with_direction(mask);
    episode_paths(i) = nansum(abs(episode_speeds)) / SpeedSamplingRate; % cm
    episode_directions(i) = sign(nanmean(speed_with_direction(mask)));
    episode_mean_speeds(i) = nanmean(abs(episode_speeds));
end

%%%%%%%%%%%%%%%%%%%%%%% COMPUTE SESSION METRICS %%%%%%%%%%%%%%%%%%%%%%%
total_session_duration = time_encoder(end) - time_encoder(1);
total_filtered_running_time = sum(episode_durations);
total_number_of_running_episodes = length(episode_starts);
total_distance_traveled = sum(episode_paths);
mean_session_speed = nanmean(speed_smth);
time_of_first_rotary_onset = pulse_times(1);
average_speed_all_episodes = mean(episode_mean_speeds, 'omitnan');

% Count total number of rewards delivered
total_number_of_rewards = length(reward_timestamps);

% Count number of episodes ending with a reward within 0.5s window
reward_window = 0.5; % seconds
episodes_with_reward = 0;
for i = 1:length(episode_ends)
    if any(reward_timestamps >= episode_ends(i) & reward_timestamps <= episode_ends(i) + reward_window)
        episodes_with_reward = episodes_with_reward + 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%% BUILD SESSION SUMMARY TABLE %%%%%%%%%%%%%%%%%%%%%%%
SessionSummary = table(...
    total_session_duration, ...
    total_filtered_running_time, ...
    total_number_of_running_episodes, ...
    total_distance_traveled, ...
    mean_session_speed, ...
    time_of_first_rotary_onset, ...
    average_speed_all_episodes, ...
    total_number_of_rewards, ...
    episodes_with_reward);


SessionSummary.Properties.VariableNames = { ...
    'Total_Session_Duration_s', ...
    'Total_Filtered_Running_Time_s', ...
    'Total_Number_of_Running_Episodes', ...
    'Total_Distance_Traveled_cm', ...
    'Mean_Session_Speed_cm_per_s', ...
    'Time_of_First_Rotary_Encoder_Onset_s', ...
    'Average_Speed_of_Running_Episodes_cm_per_s', ...
    'Total_Number_of_Rewards', ...
    'Number_of_Running_Episodes_Ending_with_Reward' ...
};


%%%%%%%%%%%%%%%%%%%%%%% LABEL DIRECTIONS AS 'CW' OR 'CCW' %%%%%%%%%%%%%%%%%%%%%%%
direction_labels = cell(size(episode_directions));
direction_labels(episode_directions == -1) = {'CW'};
direction_labels(episode_directions ==  1) = {'CCW'};
direction_labels(episode_directions ==  0) = {'Unknown'};

%%%%%%%%%%%%%%%%%%%%%%% BUILD EPISODE TABLE WITH LABELS %%%%%%%%%%%%%%%%%%%%%%%
EpisodeTable = table(...
    episode_starts, ...
    episode_durations, ...
    direction_labels, ...
    (1:length(episode_starts))', ...
    episode_paths, ...
    episode_mean_speeds);

EpisodeTable.Properties.VariableNames = { ...
    'Start_Time_s', ...
    'Duration_s', ...
    'Direction', ...
    'Episode_Number', ...
    'Total_Distance_cm', ...
    'Mean_Speed_cm_per_s' ...
};


%%%%%%%%%%%%%%%%%%%%%%% CONVERT TABLES TO CELL ARRAYS %%%%%%%%%%%%%%%%%%%%%%%
summary_cell = [SessionSummary.Properties.VariableNames; table2cell(SessionSummary)];
episode_cell = [EpisodeTable.Properties.VariableNames; table2cell(EpisodeTable)];

% Determine the maximum number of columns
max_cols = max(size(summary_cell, 2), size(episode_cell, 2));

% Pad summary_cell if needed
if size(summary_cell, 2) < max_cols
    pad_size = max_cols - size(summary_cell, 2);
    summary_cell = [summary_cell, repmat({''}, size(summary_cell, 1), pad_size)];
end

% Pad episode_cell if needed
if size(episode_cell, 2) < max_cols
    pad_size = max_cols - size(episode_cell, 2);
    episode_cell = [episode_cell, repmat({''}, size(episode_cell, 1), pad_size)];
end

% Combine into one aligned output cell
output_cell = [summary_cell; repmat({''}, 2, max_cols); episode_cell];

%%%%%%%%%%%%%%%%%%%%%%% BUILD REWARD TABLE %%%%%%%%%%%%%%%%%%%%%%%
RewardTable = table((1:length(reward_timestamps))', reward_timestamps(:));
RewardTable.Properties.VariableNames = {'Reward_Number', 'Reward_Time_s'};

%%%%%%%%%%%%%%%%%%%%%%% CONVERT REWARD TABLE TO CELL ARRAY %%%%%%%%%%%%%%%%%%%%%%%
reward_cell = [RewardTable.Properties.VariableNames; table2cell(RewardTable)];

% Determine the maximum number of columns again (in case reward_cell is larger)
max_cols = max([size(summary_cell, 2), size(episode_cell, 2), size(reward_cell, 2)]);

% Pad summary_cell if needed
if size(summary_cell, 2) < max_cols
    pad_size = max_cols - size(summary_cell, 2);
    summary_cell = [summary_cell, repmat({''}, size(summary_cell, 1), pad_size)];
end

% Pad episode_cell if needed
if size(episode_cell, 2) < max_cols
    pad_size = max_cols - size(episode_cell, 2);
    episode_cell = [episode_cell, repmat({''}, size(episode_cell, 1), pad_size)];
end

% Pad reward_cell if needed
if size(reward_cell, 2) < max_cols
    pad_size = max_cols - size(reward_cell, 2);
    reward_cell = [reward_cell, repmat({''}, size(reward_cell, 1), pad_size)];
end

% Combine all sections with spacing between them
output_cell = [summary_cell; repmat({''}, 2, max_cols); ...
               episode_cell; repmat({''}, 2, max_cols); ...
               reward_cell];

% %%%%%%%%%%%%%%%%%%%%% EXPORT TO /storage3/manus/Tables/ %%%%%%%%%%%%%%%%%%%%%%%
% ExportDir = '/storage3/manus/Tables/';
% if ~exist(ExportDir, 'dir')
%     mkdir(ExportDir);
% end
% 
% [~, trial_name] = fileparts(FileTrial);
% ExportFileName = [trial_name '.xlsx'];
% ExportPath = fullfile(ExportDir, ExportFileName);
% 
% Write summary and episode data starting at A1
% writecell([summary_cell; repmat({''}, 2, size(summary_cell,2)); episode_cell], ExportPath, 'Sheet', 1, 'Range', 'A1');
% 
% 
% Write reward headers to G5 and H5
% writecell(reward_cell(1, :), ExportPath, 'Sheet', 1, 'Range', 'G5');
% 
% Write reward data starting from G6 and H6 downwards
% writecell(reward_cell(2:end, :), ExportPath, 'Sheet', 1, 'Range', 'G6');
% 
% disp(['Exported structured summary to ', ExportPath]);

