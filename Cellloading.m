% Define the folder where the Excel files are stored
folderPath = '/storage3/manus/Tables/';  % Update this with your actual Linux path

% Get a list of all Excel files in the folder
files = dir(fullfile(folderPath, '*.xlsx'));

%%%%%%%%%%%%%%%%Single cell extraction with animal and trial info %%%%%%%%%%%%
%%%%%%%%%%%%% session duration A2 %%%%%%%%%%%%%
% Initialize a cell array to hold all A2 values (session durations)
allSessionDurations = {};

% Loop over each Excel file
for i = 1:length(files)
    % Get the full file path
    filePath = fullfile(folderPath, files(i).name);
    
    % Extract animal and trial number from the file name
    [animal, trialNumber] = extractAnimalAndTrialFromFileName(files(i).name);
    
    % Skip the file if animal or trialNumber is invalid (NaN)
    if isempty(animal) || isnan(trialNumber)
        disp(['Skipping file due to invalid animal or trial number: ', files(i).name]);
        continue;
    end
    
    % Create a key combining animal and trial number (e.g., "MP01_trial10")
    trialKey = [animal, '_trial', num2str(trialNumber)];
    
    % Get a list of all sheet names in the current Excel file
    sheetNames = sheetnames(filePath);
    
    % Loop over each sheet in the current Excel file
    for j = 1:length(sheetNames)
        % Read the specific cell (A2) from the current sheet
        sessionDuration = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'A2');
        sessionDuration = sessionDuration{1};  % Ensure it's a single value

        % Ensure sessionDuration is stored as a numeric value
        if isnumeric(sessionDuration)
            allSessionDurations{end+1} = struct('trialKey', trialKey, 'value', sessionDuration);
        elseif ischar(sessionDuration)
            numValue = str2double(sessionDuration);
            if ~isnan(numValue)
                allSessionDurations{end+1} = struct('trialKey', trialKey, 'value', numValue);
            else
                allSessionDurations{end+1} = struct('trialKey', trialKey, 'value', NaN);
            end
        else
            allSessionDurations{end+1} = struct('trialKey', trialKey, 'value', NaN);
        end
    end
end

% Display all the values in the cell array for session durations (A2)
disp('All Session Durations (from A2):');
for k = 1:length(allSessionDurations)
    fprintf('Trial %s: Session Duration: %f\n', allSessionDurations{k}.trialKey, allSessionDurations{k}.value);
end

%%%%%%%%%%%%% running time B2 %%%%%%%%%%%%%
% Initialize a cell array to hold all B2 values (running times)
allRunningTimes = {};

% Loop over each Excel file
for i = 1:length(files)
    % Get the full file path
    filePath = fullfile(folderPath, files(i).name);
    
    % Extract animal and trial number from the file name
    [animal, trialNumber] = extractAnimalAndTrialFromFileName(files(i).name);
    
    % Skip the file if animal or trialNumber is invalid (NaN)
    if isempty(animal) || isnan(trialNumber)
        disp(['Skipping file due to invalid animal or trial number: ', files(i).name]);
        continue;
    end
    
    % Create a key combining animal and trial number (e.g., "MP01_trial10")
    trialKey = [animal, '_trial', num2str(trialNumber)];
    
    % Get a list of all sheet names in the current Excel file
    sheetNames = sheetnames(filePath);
    
    % Loop over each sheet in the current Excel file
    for j = 1:length(sheetNames)
        % Read the specific cell (B2) from the current sheet
        runningTime = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'B2');
        runningTime = runningTime{1};  % Ensure it's a single value

        % Check if the runningTime is numeric
        if isnumeric(runningTime)
            allRunningTimes{end+1} = struct('trialKey', trialKey, 'value', runningTime);
        elseif ischar(runningTime)
            numValue = str2double(runningTime);
            if ~isnan(numValue)
                allRunningTimes{end+1} = struct('trialKey', trialKey, 'value', numValue);
            else
                allRunningTimes{end+1} = struct('trialKey', trialKey, 'value', runningTime);
            end
        else
            allRunningTimes{end+1} = struct('trialKey', trialKey, 'value', runningTime);
        end
    end
end

% Display the collected B2 values (running times) using fprintf
disp('All Total Running Times (from B2):');
for k = 1:length(allRunningTimes)
    fprintf('Trial %s: Running Time: %f\n', allRunningTimes{k}.trialKey, allRunningTimes{k}.value);
end

%%%%%%%%%%%%% total number of running episodes C2 %%%%%%%%%%%%%
% Initialize a cell array to hold all C2 values (number of running episodes)
allRunningEpisodes = {};

% Loop over each Excel file
for i = 1:length(files)
    % Get the full file path
    filePath = fullfile(folderPath, files(i).name);
    
    % Extract animal and trial number from the file name
    [animal, trialNumber] = extractAnimalAndTrialFromFileName(files(i).name);
    
    % Skip the file if animal or trialNumber is invalid (NaN)
    if isempty(animal) || isnan(trialNumber)
        disp(['Skipping file due to invalid animal or trial number: ', files(i).name]);
        continue;
    end
    
    % Create a key combining animal and trial number (e.g., "MP01_trial10")
    trialKey = [animal, '_trial', num2str(trialNumber)];
    
    % Get a list of all sheet names in the current Excel file
    sheetNames = sheetnames(filePath);
    
    % Loop over each sheet in the current Excel file
    for j = 1:length(sheetNames)
        % Read the specific cell (C2) from the current sheet
        runningEpisodes = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'C2');
        runningEpisodes = runningEpisodes{1};  % Ensure it's a single value

        % Ensure runningEpisodes is stored as a numeric value
        if isnumeric(runningEpisodes)
            allRunningEpisodes{end+1} = struct('trialKey', trialKey, 'value', runningEpisodes);
        elseif ischar(runningEpisodes)
            numValue = str2double(runningEpisodes);
            if ~isnan(numValue)
                allRunningEpisodes{end+1} = struct('trialKey', trialKey, 'value', numValue);
            else
                allRunningEpisodes{end+1} = struct('trialKey', trialKey, 'value', NaN);
            end
        else
            allRunningEpisodes{end+1} = struct('trialKey', trialKey, 'value', NaN);
        end
    end
end

% Display the collected C2 values (total number of running episodes) using fprintf
disp('All Total Number of Running Episodes (from C2):');
for k = 1:length(allRunningEpisodes)
    fprintf('Trial %s: Running Episode Count: %f\n', allRunningEpisodes{k}.trialKey, allRunningEpisodes{k}.value);
end

%%%%%%%%%%%%% total distance traveled D2 %%%%%%%%%%%%%
% Initialize a cell array to hold all D2 values (distance traveled)
allDistanceTraveled = {};

% Loop over each Excel file
for i = 1:length(files)
    % Get the full file path
    filePath = fullfile(folderPath, files(i).name);
    
    % Extract animal and trial number from the file name
    [animal, trialNumber] = extractAnimalAndTrialFromFileName(files(i).name);
    
    % Skip the file if animal or trialNumber is invalid (NaN)
    if isempty(animal) || isnan(trialNumber)
        disp(['Skipping file due to invalid animal or trial number: ', files(i).name]);
        continue;
    end
    
    % Create a key combining animal and trial number (e.g., "MP01_trial10")
    trialKey = [animal, '_trial', num2str(trialNumber)];
    
    % Get a list of all sheet names in the current Excel file
    sheetNames = sheetnames(filePath);
    
    % Loop over each sheet in the current Excel file
    for j = 1:length(sheetNames)
        % Read the specific cell (D2) from the current sheet
        distanceTraveled = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'D2');
        distanceTraveled = distanceTraveled{1};  % Ensure it's a single value

        % Ensure distanceTraveled is stored as a numeric value
        if isnumeric(distanceTraveled)
            allDistanceTraveled{end+1} = struct('trialKey', trialKey, 'value', distanceTraveled);
        elseif ischar(distanceTraveled)
            numValue = str2double(distanceTraveled);
            if ~isnan(numValue)
                allDistanceTraveled{end+1} = struct('trialKey', trialKey, 'value', numValue);
            else
                allDistanceTraveled{end+1} = struct('trialKey', trialKey, 'value', NaN);
            end
        else
            allDistanceTraveled{end+1} = struct('trialKey', trialKey, 'value', NaN);
        end
    end
end

% Display the collected D2 values (total distance traveled) using fprintf
disp('All Total Distance Traveled (from D2):');
for k = 1:length(allDistanceTraveled)
    fprintf('Trial %s: Distance Traveled: %f\n', allDistanceTraveled{k}.trialKey, allDistanceTraveled{k}.value);
end

%%%%%%%%%%%%% mean session speed E2 %%%%%%%%%%%%%
% Initialize a cell array to hold all E2 values (mean session speed)
allMeanSessionSpeed = {};

% Loop over each Excel file
for i = 1:length(files)
    % Get the full file path
    filePath = fullfile(folderPath, files(i).name);
    
    % Extract animal and trial number from the file name
    [animal, trialNumber] = extractAnimalAndTrialFromFileName(files(i).name);
    
    % Skip the file if animal or trialNumber is invalid (NaN)
    if isempty(animal) || isnan(trialNumber)
        disp(['Skipping file due to invalid animal or trial number: ', files(i).name]);
        continue;
    end
    
    % Create a key combining animal and trial number (e.g., "MP01_trial10")
    trialKey = [animal, '_trial', num2str(trialNumber)];
    
    % Get a list of all sheet names in the current Excel file
    sheetNames = sheetnames(filePath);
    
    % Loop over each sheet in the current Excel file
    for j = 1:length(sheetNames)
        % Read the specific cell (E2) from the current sheet
        meanSessionSpeed = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'E2');
        meanSessionSpeed = meanSessionSpeed{1};  % Ensure it's a single value

        % Ensure meanSessionSpeed is stored as a numeric value
        if isnumeric(meanSessionSpeed)
            allMeanSessionSpeed{end+1} = struct('trialKey', trialKey, 'value', meanSessionSpeed);
        elseif ischar(meanSessionSpeed)
            numValue = str2double(meanSessionSpeed);
            if ~isnan(numValue)
                allMeanSessionSpeed{end+1} = struct('trialKey', trialKey, 'value', numValue);
            else
                allMeanSessionSpeed{end+1} = struct('trialKey', trialKey, 'value', NaN);
            end
        else
            allMeanSessionSpeed{end+1} = struct('trialKey', trialKey, 'value', NaN);
        end
    end
end

% Display the collected E2 values (mean session speed) using fprintf
disp('All Mean Session Speed (from E2):');
for k = 1:length(allMeanSessionSpeed)
    fprintf('Trial %s: Mean Session Speed: %f\n', allMeanSessionSpeed{k}.trialKey, allMeanSessionSpeed{k}.value);
end

%%%%%%%%%%%%% time of first rotary encoder onset F2 %%%%%%%%%%%%%
% Initialize a cell array to hold all F2 values (first rotary encoder onset time)
allRotaryEncoderOnset = {};

% Loop over each Excel file
for i = 1:length(files)
    % Get the full file path
    filePath = fullfile(folderPath, files(i).name);
    
    % Extract animal and trial number from the file name
    [animal, trialNumber] = extractAnimalAndTrialFromFileName(files(i).name);
    
    % Skip the file if animal or trialNumber is invalid (NaN)
    if isempty(animal) || isnan(trialNumber)
        disp(['Skipping file due to invalid animal or trial number: ', files(i).name]);
        continue;
    end
    
    % Create a key combining animal and trial number (e.g., "MP01_trial10")
    trialKey = [animal, '_trial', num2str(trialNumber)];
    
    % Get a list of all sheet names in the current Excel file
    sheetNames = sheetnames(filePath);
    
    % Loop over each sheet in the current Excel file
    for j = 1:length(sheetNames)
        % Read the specific cell (F2) from the current sheet
        rotaryEncoderOnset = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'F2');
        rotaryEncoderOnset = rotaryEncoderOnset{1};  % Ensure it's a single value

        % Ensure rotaryEncoderOnset is stored as a numeric value
        if isnumeric(rotaryEncoderOnset)
            allRotaryEncoderOnset{end+1} = struct('trialKey', trialKey, 'value', rotaryEncoderOnset);
        elseif ischar(rotaryEncoderOnset)
            numValue = str2double(rotaryEncoderOnset);
            if ~isnan(numValue)
                allRotaryEncoderOnset{end+1} = struct('trialKey', trialKey, 'value', numValue);
            else
                allRotaryEncoderOnset{end+1} = struct('trialKey', trialKey, 'value', NaN);
            end
        else
            allRotaryEncoderOnset{end+1} = struct('trialKey', trialKey, 'value', NaN);
        end
    end
end

% Display the collected F2 values (time of first rotary encoder onset) using fprintf
disp('All Time of First Rotary Encoder Onset (from F2):');
for k = 1:length(allRotaryEncoderOnset)
    fprintf('Trial %s: Rotary Encoder Onset: %f\n', allRotaryEncoderOnset{k}.trialKey, allRotaryEncoderOnset{k}.value);
end

%%%%%%%%%%%%% average speed of running episodes G2 %%%%%%%%%%%%%
% Initialize a cell array to hold all G2 values (average speed of running episodes)
allAvgRunningSpeed = {};

% Loop over each Excel file
for i = 1:length(files)
    % Get the full file path
    filePath = fullfile(folderPath, files(i).name);
    
    % Extract animal and trial number from the file name
    [animal, trialNumber] = extractAnimalAndTrialFromFileName(files(i).name);
    
    % Skip the file if animal or trialNumber is invalid (NaN)
    if isempty(animal) || isnan(trialNumber)
        disp(['Skipping file due to invalid animal or trial number: ', files(i).name]);
        continue;
    end
    
    % Create a key combining animal and trial number (e.g., "MP01_trial10")
    trialKey = [animal, '_trial', num2str(trialNumber)];
    
    % Get a list of all sheet names in the current Excel file
    sheetNames = sheetnames(filePath);
    
    % Loop over each sheet in the current Excel file
    for j = 1:length(sheetNames)
        % Read the specific cell (G2) from the current sheet
        avgRunningSpeed = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'G2');
        avgRunningSpeed = avgRunningSpeed{1};  % Ensure it's a single value

        % Ensure avgRunningSpeed is stored as a numeric value
        if isnumeric(avgRunningSpeed)
            allAvgRunningSpeed{end+1} = struct('trialKey', trialKey, 'value', avgRunningSpeed);
        elseif ischar(avgRunningSpeed)
            numValue = str2double(avgRunningSpeed);
            if ~isnan(numValue)
                allAvgRunningSpeed{end+1} = struct('trialKey', trialKey, 'value', numValue);
            else
                allAvgRunningSpeed{end+1} = struct('trialKey', trialKey, 'value', NaN);
            end
        else
            allAvgRunningSpeed{end+1} = struct('trialKey', trialKey, 'value', NaN);
        end
    end
end

% Display the collected G2 values (average speed of running episodes) using fprintf
disp('All Average Speed of Running Episodes (from G2):');
for k = 1:length(allAvgRunningSpeed)
    fprintf('Trial %s: Avg Running Speed: %f\n', allAvgRunningSpeed{k}.trialKey, allAvgRunningSpeed{k}.value);
end

%%%%%%%%%%%%% total number of rewards H2 %%%%%%%%%%%%%
% Initialize a cell array to hold all H2 values (total number of rewards)
allRewards = {};

% Loop over each Excel file
for i = 1:length(files)
    % Get the full file path
    filePath = fullfile(folderPath, files(i).name);
    
    % Extract animal and trial number from the file name
    [animal, trialNumber] = extractAnimalAndTrialFromFileName(files(i).name);
    
    % Skip the file if animal or trialNumber is invalid (NaN)
    if isempty(animal) || isnan(trialNumber)
        disp(['Skipping file due to invalid animal or trial number: ', files(i).name]);
        continue;
    end
    
    % Create a key combining animal and trial number (e.g., "MP01_trial10")
    trialKey = [animal, '_trial', num2str(trialNumber)];
    
    % Get a list of all sheet names in the current Excel file
    sheetNames = sheetnames(filePath);
    
    % Loop over each sheet in the current Excel file
    for j = 1:length(sheetNames)
        % Read the specific cell (H2) from the current sheet
        rewards = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'H2');
        rewards = rewards{1};  % Ensure it's a single value

        % Ensure rewards is stored as a numeric value
        if isnumeric(rewards)
            allRewards{end+1} = struct('trialKey', trialKey, 'value', rewards);
        elseif ischar(rewards)
            numValue = str2double(rewards);
            if ~isnan(numValue)
                allRewards{end+1} = struct('trialKey', trialKey, 'value', numValue);
            else
                allRewards{end+1} = struct('trialKey', trialKey, 'value', NaN);
            end
        else
            allRewards{end+1} = struct('trialKey', trialKey, 'value', NaN);
        end
    end
end

% Display the collected H2 values (total number of rewards) using fprintf
disp('All Total Number of Rewards (from H2):');
for k = 1:length(allRewards)
    fprintf('Trial %s: Reward Count: %f\n', allRewards{k}.trialKey, allRewards{k}.value);
end

%%%%%%%%%%%%% body weight J2 %%%%%%%%%%%%%
% Initialize a cell array to hold all J2 values (body weights)
allBodyWeights = {};

% Loop over each Excel file
for i = 1:length(files)
    % Get the full file path
    filePath = fullfile(folderPath, files(i).name);
    
    % Extract animal and trial number from the file name
    [animal, trialNumber] = extractAnimalAndTrialFromFileName(files(i).name);
    
    % Skip the file if animal or trialNumber is invalid (NaN)
    if isempty(animal) || isnan(trialNumber)
        disp(['Skipping file due to invalid animal or trial number: ', files(i).name]);
        continue;
    end
    
    % Create a key combining animal and trial number (e.g., "MP01_trial10")
    trialKey = [animal, '_trial', num2str(trialNumber)];
    
    % Get a list of all sheet names in the current Excel file
    sheetNames = sheetnames(filePath);
    
    % Loop over each sheet in the current Excel file
    for j = 1:length(sheetNames)
        % Read the specific cell (J2) from the current sheet
        bodyWeight = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'J2');
        bodyWeight = bodyWeight{1};  % Ensure it's a single value

        % Ensure bodyWeight is stored as a numeric value
        if isnumeric(bodyWeight)
            allBodyWeights{end+1} = struct('trialKey', trialKey, 'value', bodyWeight);
        elseif ischar(bodyWeight)
            numValue = str2double(bodyWeight);
            if ~isnan(numValue)
                allBodyWeights{end+1} = struct('trialKey', trialKey, 'value', numValue);
            else
                allBodyWeights{end+1} = struct('trialKey', trialKey, 'value', NaN);
            end
        else
            allBodyWeights{end+1} = struct('trialKey', trialKey, 'value', NaN);
        end
    end
end

% Display the collected J2 values (body weights) using fprintf
disp('All Body Weights (from J2):');
for k = 1:length(allBodyWeights)
    fprintf('Trial %s: Body Weight: %f\n', allBodyWeights{k}.trialKey, allBodyWeights{k}.value);
end


%%%%%%%%%%%%%%%% multiple cell extraction %%%%%%%%%%%5
% Initialize a containers.Map to hold the data grouped by animal and trial
% Initialize a containers.Map to hold the data grouped by animal and trial
allEpisodeDataByTrial = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Loop over each Excel file
for i = 1:length(files)
    % Get the full file path
    filePath = fullfile(folderPath, files(i).name);
    
    % Extract animal and trial number from the file name
    [animal, trialNumber] = extractAnimalAndTrialFromFileName(files(i).name);
    
    % Skip the file if animal or trialNumber is invalid (NaN)
    if isempty(animal) || isnan(trialNumber)
        disp(['Skipping file due to invalid animal or trial number: ', files(i).name]);
        continue;
    end
    
    % Create a key combining animal and trial number (e.g., "MP01_trial10")
    trialKey = [animal, '_trial', num2str(trialNumber)];
    
    % Get a list of all sheet names in the current Excel file
    sheetNames = sheetnames(filePath);
    
    % Loop over each sheet in the current Excel file
    for j = 1:length(sheetNames)
        % Read the range of cells for the parameters starting from row 6
        episodeStartTimes = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'A6:A100'); % Start times from column A
        episodeLengths = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'B6:B100'); % Episode lengths from column B
        directions = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'C6:C100'); % Directions from column C
        episodeNumbers = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'D6:D100'); % Running episode numbers from column D
        distances = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'E6:E100'); % Distances from column E
        meanSpeeds = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'F6:F100'); % Mean speeds from column F
        
        % Ensure the data read is not empty and is in the expected format
        if isempty(episodeStartTimes)
            disp(['No data found in A6:A100 for ', files(i).name, ', Sheet: ', sheetNames{j}]);
            continue;
        end
        
        % Initialize the trialData structure for the current trial if not already initialized
        if ~isKey(allEpisodeDataByTrial, trialKey)
            % Initialize the structure if it does not exist yet
            allEpisodeDataByTrial(trialKey) = struct( ...
    'startTimes', {{}}, ...
    'episodeLengths', {{}}, ...
    'directions', {{}}, ...
    'episodeNumbers', {{}}, ...
    'distances', {{}}, ...
    'meanSpeeds', {{}} );
        end
        
        % Access the current trial data
        trialData = allEpisodeDataByTrial(trialKey);
        
        % Ensure all fields in trialData are initialized (if they are not already initialized)
        if ~isfield(trialData, 'startTimes')
            trialData.startTimes = {}; % Initialize if empty
        end
        if ~isfield(trialData, 'episodeLengths')
            trialData.episodeLengths = {}; % Initialize if empty
        end
        if ~isfield(trialData, 'directions')
            trialData.directions = {}; % Initialize if empty
        end
        if ~isfield(trialData, 'episodeNumbers')
            trialData.episodeNumbers = {}; % Initialize if empty
        end
        if ~isfield(trialData, 'distances')
            trialData.distances = {}; % Initialize if empty
        end
        if ~isfield(trialData, 'meanSpeeds')
            trialData.meanSpeeds = {}; % Initialize if empty
        end
        
        % Loop through the read values and store them in the trial's list
        for k = 1:length(episodeStartTimes)
            % Extract the corresponding values for this episode
            startTime = episodeStartTimes{k};
            episodeLength = episodeLengths{k};
            direction = directions{k};
            episodeNumber = episodeNumbers{k};
            distance = distances{k};
            meanSpeed = meanSpeeds{k};
            
            % Ensure that data is not empty or NaN before appending to the trialData
            if ~isempty(startTime) && ~isnan(startTime) && ...
               ~isempty(episodeLength) && ~isnan(episodeLength) && ...
               ~isempty(direction) && ~isnan(episodeNumber) && ...
               ~isempty(distance) && ~isnan(meanSpeed)
                % Append the data to the trial's structure
                trialData.startTimes{end+1} = startTime;
                trialData.episodeLengths{end+1} = episodeLength;
                trialData.directions{end+1} = direction;
                trialData.episodeNumbers{end+1} = episodeNumber;
                trialData.distances{end+1} = distance;
                trialData.meanSpeeds{end+1} = meanSpeed;
                
                % Update the map entry
                allEpisodeDataByTrial(trialKey) = trialData;
            else
                % If any value is invalid, break and skip to the next episode
                disp(['Skipping invalid episode data for trial ', trialKey, ', Episode: ', num2str(k)]);
                break;
            end
        end
    end
end

% Display the collected episode data, grouped by animal and trial
disp('All Episode Data Grouped by Animal and Trial:');
trialKeys = keys(allEpisodeDataByTrial);
for i = 1:length(trialKeys)
    trial = trialKeys{i};
    disp(['Trial ', trial, ' Episode Data:']);
    
    episodeData = allEpisodeDataByTrial(trial);
    for k = 1:length(episodeData.startTimes)
        fprintf('  Episode %d:\n', k);
        fprintf('    Start Time: %f\n', episodeData.startTimes{k});
        fprintf('    Length: %f\n', episodeData.episodeLengths{k});
        fprintf('    Direction: %s\n', episodeData.directions{k});
        fprintf('    Episode Number: %d\n', episodeData.episodeNumbers{k});
        fprintf('    Distance: %f\n', episodeData.distances{k});
        fprintf('    Mean Speed: %f\n', episodeData.meanSpeeds{k});
    end
end

%%%%%%%%%%%%% rewards and reward times G6 and H6 %%%%%%%%%%%%%
%%%%%%%%%%%%% rewards and reward times G6 and H6 %%%%%%%%%%%%%
% Initialize a containers.Map to hold the reward data grouped by animal and trial
allRewardsDataByTrial = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Loop over each Excel file
for i = 1:length(files)
    % Get the full file path
    filePath = fullfile(folderPath, files(i).name);
    
    % Extract animal and trial number from the file name
    [animal, trialNumber] = extractAnimalAndTrialFromFileName(files(i).name);
    
    % Skip the file if animal or trialNumber is invalid (NaN)
    if isempty(animal) || isnan(trialNumber)
        disp(['Skipping file due to invalid animal or trial number: ', files(i).name]);
        continue;
    end
    
    % Create a key combining animal and trial number (e.g., "MP01_trial10")
    trialKey = [animal, '_trial', num2str(trialNumber)];
    
    % Get a list of all sheet names in the current Excel file
    sheetNames = sheetnames(filePath);
    
    % Loop over each sheet in the current Excel file
    for j = 1:length(sheetNames)
        % Read the range of cells for the rewards (column G) and reward times (column H)
        rewardCounts = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'G6:G100'); % Rewards count from column G
        rewardTimes = readcell(filePath, 'Sheet', sheetNames{j}, 'Range', 'H6:H100'); % Reward times from column H
        
        % Ensure the data read is not empty and is in the expected format
        if isempty(rewardCounts) || isempty(rewardTimes)
            disp(['No data found for rewards or reward times in file ', files(i).name, ', Sheet: ', sheetNames{j}]);
            continue;
        end
        
        % Initialize the trialData structure for the current trial if not already initialized
        if ~isKey(allRewardsDataByTrial, trialKey)
            % Initialize the structure if it does not exist yet
            allRewardsDataByTrial(trialKey) = struct( ...
    'rewardCounts', {{}}, ...
    'rewardTimes', {{}});
        end
        
        % Access the current trial data
        trialData = allRewardsDataByTrial(trialKey);
        
        % Ensure all fields in trialData are initialized (if they are not already initialized)
        if ~isfield(trialData, 'rewardCounts')
            trialData.rewardCounts = {}; % Initialize if empty
        end
        if ~isfield(trialData, 'rewardTimes')
            trialData.rewardTimes = {}; % Initialize if empty
        end
        
        % Loop through the read values and store them in the trial's list
        for k = 1:length(rewardCounts)
            % Extract the corresponding values for this reward
            rewardCount = rewardCounts{k};
            rewardTime = rewardTimes{k};
            
            % Check if the reward count is missing or empty, break if so
            if ismissing(rewardCount) || isempty(rewardCount)
                break; % No more rewards after this point
            end
            
            % Ensure that both reward count and reward time are valid before appending
            if ~ismissing(rewardCount) && ~isempty(rewardCount) && ...
               ~ismissing(rewardTime) && ~isempty(rewardTime)
                % Append the reward data to the trial's structure
                trialData.rewardCounts{end+1} = rewardCount;
                trialData.rewardTimes{end+1} = rewardTime;
                
                % Update the map entry
                allRewardsDataByTrial(trialKey) = trialData;
            else
                % If any value is invalid, break and skip to the next reward
                disp(['Skipping invalid reward data for trial ', trialKey, ', Reward: ', num2str(k)]);
                break;
            end
        end
    end
end

% Display the collected reward data, grouped by animal and trial
disp('All Reward Data Grouped by Animal and Trial:');
trialKeys = keys(allRewardsDataByTrial);
for i = 1:length(trialKeys)
    trial = trialKeys{i};
    disp(['Trial ', trial, ' Reward Data:']);
    
    rewardData = allRewardsDataByTrial(trial);
    for k = 1:length(rewardData.rewardCounts)
        fprintf('  Reward %d:\n', k);
        fprintf('    Reward Count: %d\n', rewardData.rewardCounts{k});
        fprintf('    Reward Time: %f\n', rewardData.rewardTimes{k});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULATIONS%%%%%%%%%%%%%%%%%%%%



% Define starting weights for each animal
startingWeights = containers.Map({'MP01', 'MP02'}, {610, 470});  % MP01 starts at 610g, MP02 starts at 470g

% Initialize a cell array to hold the body weight percentage data
allBodyWeightPercentages = {};

% Loop over each file (assuming you have a similar loop for extracting trial information)
for i = 1:length(files)
    % Get the full file path
    filePath = fullfile(folderPath, files(i).name);
    
    % Extract animal and trial number from the file name
    [animal, trialNumber] = extractAnimalAndTrialFromFileName(files(i).name);
    
    % Skip the file if animal or trialNumber is invalid (NaN)
    if isempty(animal) || isnan(trialNumber)
        disp(['Skipping file due to invalid animal or trial number: ', files(i).name]);
        continue;
    end
    
    % Create a key combining animal and trial number (e.g., "MP01_trial10")
    trialKey = [animal, '_trial', num2str(trialNumber)];
    
    % Get the body weight (J2) for this trial
    bodyWeight = readcell(filePath, 'Sheet', sheetNames{1}, 'Range', 'J2');  % Assuming it’s from the first sheet
    bodyWeight = bodyWeight{1};  % Ensure it's a single value
    
    % Retrieve the starting weight from the map for the current animal
    if isKey(startingWeights, animal)
        startingWeight = startingWeights(animal);
    else
        disp(['Starting weight for ', animal, ' is not available. Skipping trial.']);
        continue;
    end
    
    % Calculate the percentage of the starting weight
    weightPercentage = (bodyWeight / startingWeight) * 100;
    
    % Store the result in the array
    allBodyWeightPercentages{end+1} = struct('trialKey', trialKey, 'bodyWeight', bodyWeight, 'percentage', weightPercentage);
end

% Display the collected body weight percentages
disp('All Body Weight Percentages:');
for k = 1:length(allBodyWeightPercentages)
    fprintf('Trial %s: Body Weight: %f g, Percentage of Starting Weight: %f%%\n', ...
            allBodyWeightPercentages{k}.trialKey, allBodyWeightPercentages{k}.bodyWeight, ...
            allBodyWeightPercentages{k}.percentage);
end

