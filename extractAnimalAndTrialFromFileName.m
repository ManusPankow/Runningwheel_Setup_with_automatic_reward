% Function to extract the animal and trial number from the file name
function [animal, trialNumber] = extractAnimalAndTrialFromFileName(fileName)
    % Use regular expression to extract the animal identifier (MP01 or MP02) and trial number (e.g., trial10)
    animalPattern = '^(MP01|MP02)';  % Look for MP01 or MP02 at the beginning
    trialPattern = 'trial(\d+)';     % Look for trial number in the format 'trial10', 'trial1', etc.

    animalMatch = regexp(fileName, animalPattern, 'match');
    trialMatch = regexp(fileName, trialPattern, 'tokens');

    % If animal is found, extract it; otherwise, return empty
    if ~isempty(animalMatch)
        animal = animalMatch{1};
    else
        animal = '';  % Return empty if no animal is found
    end

    % If trial number is found, convert it to numeric; otherwise, return NaN
    if ~isempty(trialMatch)
        trialNumber = str2double(trialMatch{1}{1});
    else
        trialNumber = NaN;  % Return NaN if no trial number is found
    end
end