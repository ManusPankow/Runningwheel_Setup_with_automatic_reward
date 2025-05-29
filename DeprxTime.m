% Assuming allBodyWeights and allRunningTimes are available

% Create dictionaries for session durations and running times
session_durations_dict = containers.Map();
running_times_dict = containers.Map();

for i = 1:length(allBodyWeights)
    session_durations_dict(allBodyWeights{i}.trialKey) = allBodyWeights{i}.value;
end
for i = 1:length(allRunningTimes)
    running_times_dict(allRunningTimes{i}.trialKey) = allRunningTimes{i}.value;
end

% Get percent running time
percent_running_time = [];
trial_keys = {};
for trial_key = keys(session_durations_dict)
    trial_key = trial_key{1};
    session_duration = session_durations_dict(trial_key);
    if session_duration > 0
        percent_running_time = [percent_running_time; session_duration];
        trial_keys = [trial_keys; trial_key];
    end
end

% Separate and sort by animal
animal_1_trials = trial_keys(startsWith(trial_keys, 'MP01'));
animal_2_trials = trial_keys(startsWith(trial_keys, 'MP02'));

trial_numbers_animal_1 = cellfun(@(x) str2double(regexp(x, 'trial(\d+)', 'tokens', 'once')), animal_1_trials);
trial_numbers_animal_2 = cellfun(@(x) str2double(regexp(x, 'trial(\d+)', 'tokens', 'once')), animal_2_trials);

[trial_x_animal_1, idx1] = sort(trial_numbers_animal_1);
[trial_x_animal_2, idx2] = sort(trial_numbers_animal_2);
animal_1_trials_sorted = animal_1_trials(idx1);
animal_2_trials_sorted = animal_2_trials(idx2);

percent_animal_1 = percent_running_time(ismember(trial_keys, animal_1_trials_sorted));
percent_animal_2 = percent_running_time(ismember(trial_keys, animal_2_trials_sorted));
percent_animal_1 = percent_animal_1(idx1) / 610 * 100;
percent_animal_2 = percent_animal_2(idx2) / 470 * 100;

% Running time data
total_running_time = [];
running_keys = {};
for trial_key = keys(running_times_dict)
    trial_key = trial_key{1};
    val = running_times_dict(trial_key);
    if val > 0
        total_running_time = [total_running_time; val];
        running_keys = [running_keys; trial_key];
    end
end

animal_1_rt_trials = running_keys(startsWith(running_keys, 'MP01'));
animal_2_rt_trials = running_keys(startsWith(running_keys, 'MP02'));

session_numbers_animal_1 = cellfun(@(x) str2double(regexp(x, 'trial(\d+)', 'tokens', 'once')), animal_1_rt_trials);
session_numbers_animal_2 = cellfun(@(x) str2double(regexp(x, 'trial(\d+)', 'tokens', 'once')), animal_2_rt_trials);

[session_x_animal_1, idx_rt1] = sort(session_numbers_animal_1);
[session_x_animal_2, idx_rt2] = sort(session_numbers_animal_2);

animal_1_rt_sorted = animal_1_rt_trials(idx_rt1);
animal_2_rt_sorted = animal_2_rt_trials(idx_rt2);

running_time_animal_1 = total_running_time(ismember(running_keys, animal_1_rt_sorted));
running_time_animal_2 = total_running_time(ismember(running_keys, animal_2_rt_sorted));
running_time_animal_1 = running_time_animal_1(idx_rt1);
running_time_animal_2 = running_time_animal_2(idx_rt2);

%% Plotting
figure('Position', [100, 100, 1600, 600]);

%% Subplot 1: Deprivation
subplot(1, 2, 1); hold on;

% Axis limits
all_percent = [percent_animal_1; percent_animal_2];
ylim_vals = [floor(min(all_percent) - 2), ceil(max(all_percent) + 2)];
xmax1 = max([trial_x_animal_1; trial_x_animal_2]);

% Shading (grayscale only)
fill([0.5 4.5 4.5 0.5], [ylim_vals(1) ylim_vals(1) ylim_vals(2) ylim_vals(2)], [0.85 0.85 0.85], 'EdgeColor', 'none'); % Light gray
fill([4.5 8.5 8.5 4.5], [ylim_vals(1) ylim_vals(1) ylim_vals(2) ylim_vals(2)], [1 1 1], 'EdgeColor', 'none');           % White
fill([8.5 10.5 10.5 8.5], [ylim_vals(1) ylim_vals(1) ylim_vals(2) ylim_vals(2)], [0.85 0.85 0.85], 'EdgeColor', 'none');% Light gray
fill([10.5 xmax1+0.5 xmax1+0.5 10.5], [ylim_vals(1) ylim_vals(1) ylim_vals(2) ylim_vals(2)], [0.3 0.3 0.3], 'EdgeColor', 'none'); % Dark gray

% Plot with colorful lines
h1 = plot(trial_x_animal_1, percent_animal_1, '-o', 'LineWidth', 3, 'Color', 'b', 'MarkerSize', 12);  % Blue
h2 = plot(trial_x_animal_2, percent_animal_2, '-^', 'LineWidth', 3, 'Color', [1, 0.443, 0.122], 'MarkerSize', 12);  % Orange

xlabel('Sessions', 'FontSize', 20);
ylabel('Deprivation depth (% of initial BW)', 'FontSize', 20);
set(gca, 'FontSize', 20);
ylim(ylim_vals);
legend([h1 h2], {'MP01', 'MP02'}, 'Location', 'northwest');

%% Subplot 2: Running Time
subplot(1, 2, 2); hold on;

% Axis limits
all_rt = [running_time_animal_1; running_time_animal_2];
ylim_vals_rt = [floor(min(all_rt) - 10), ceil(max(all_rt) + 10)];
xmax2 = max([session_x_animal_1; session_x_animal_2]);

% Shading
fill([0.5 4.5 4.5 0.5], [ylim_vals_rt(1) ylim_vals_rt(1) ylim_vals_rt(2) ylim_vals_rt(2)], [0.85 0.85 0.85], 'EdgeColor', 'none');     
fill([4.5 8.5 8.5 4.5], [ylim_vals_rt(1) ylim_vals_rt(1) ylim_vals_rt(2) ylim_vals_rt(2)], [1 1 1], 'EdgeColor', 'none'); 
fill([8.5 10.5 10.5 8.5], [ylim_vals_rt(1) ylim_vals_rt(1) ylim_vals_rt(2) ylim_vals_rt(2)], [0.85 0.85 0.85], 'EdgeColor', 'none');    
fill([10.5 xmax2+0.5 xmax2+0.5 10.5], [ylim_vals_rt(1) ylim_vals_rt(1) ylim_vals_rt(2) ylim_vals_rt(2)], [0.3 0.3 0.3], 'EdgeColor', 'none'); 

% Plot with colorful lines
h3 = plot(session_x_animal_1, running_time_animal_1, '-o', 'LineWidth', 3, 'Color', 'b', 'MarkerSize', 12);  % Blue
h4 = plot(session_x_animal_2, running_time_animal_2, '-^', 'LineWidth', 3, 'Color', [1, 0.443, 0.122], 'MarkerSize', 12);  % Orange

xlabel('Sessions', 'FontSize', 20);
ylabel('Total running time (min)', 'FontSize', 20);
set(gca, 'FontSize', 20);
ylim(ylim_vals_rt);
legend([h3 h4], {'MP01', 'MP02'}, 'Location', 'northwest');

% Export to PDF
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', ...
         'PaperUnits', 'Inches', ...
         'PaperSize', [pos(3), pos(4)]);
print(gcf, 'results_subplot_final_colorlines_grayshading', '-dpdf', '-r0');
