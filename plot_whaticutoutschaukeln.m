figure('Name', 'Speed and Running Periods Analysis', 'NumberTitle', 'off');

%% Define time window
time_window = [0, 500];

%%%%%%%%%%%%%%%%%%%%%% Subplot 1: Speed with inactive periods %%%%%%%%%%%%%%%%%%%%%%%
% subplot(3,1,1);
% 
% plot(time_speed2, speed_interp, 'k-'); 
% hold on;
% plot(time_speed2, speed_smth, 'b', 'LineWidth', 1.5);
% 
% inactive_periods = [inactivity_start(:), inactivity_end(:)];
% for i = 1:size(inactive_periods, 
%     fill([inactive_periods(i,1), inactive_periods(i,2), inactive_periods(i,2), inactive_periods(i,1)], ...
%          [min(speed_interp), min(speed_interp), 120, 120], ...
%          [0.7 0.7 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% end
% 
% threshold_value = 37.72;
% yline(threshold_value, '--k', 'LineWidth', 1.5);
% 
% xlims = xlim;
% x_offset = xlims(1) + 10;  % Adjust as needed
% text(x_offset, threshold_value + 2, sprintf('%.2f cm/s', threshold_value), ...
%     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
%     'FontSize', 20, 'FontWeight', 'bold');
% 
% xlabel('Time (s)', 'FontSize', 20);
% ylabel('Speed (cm/s)', 'FontSize', 20);
% grid on;
% xlim(time_window);
% set(gca, 'FontSize', 20);
% hold off;

%%%%%%%%%%%%%%%%%%%%%% Subplot 2: Running Direction (Signed Speed) %%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,1);
hold on;

h_unsigned = plot(time_speed2, speed_smth, '.-', 'Color', [0.7 0.7 0.7], 'MarkerSize', 1);
h_signed = plot(time_speed2, speed_with_direction, '.-', 'Color', 'k', 'MarkerSize', 1, 'LineWidth', 1.5);

yline(0, '--k', 'LineWidth', 1.5, 'FontSize', 20);

xlims = xlim;
x_pos = xlims(1) + 10;  % Shifted to the right
text(x_pos,  60,  'CCW', 'HorizontalAlignment','left', 'FontWeight','bold', 'FontSize', 20);
text(x_pos, -40,  'CW',  'HorizontalAlignment','left', 'FontWeight','bold', 'FontSize', 20);

ylim([-75, 120]);

xlabel('Time (s)', 'FontSize', 20);
ylabel('Speed (cm/s)', 'FontSize', 20);
grid on;
set(gca, 'FontSize', 20);
box on
hold off;
xlim([0,480]);


%%%%%%%%%%%%%%%%%%%%%% Subplot 3: Filtered Running Periods %%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2);
plot(time_speed2, speed_smth, 'k', 'LineWidth', 1.5);
hold on;

for i = 1:size(filtered_running_periods, 1)
    period_start = filtered_running_periods(i, 1);
    period_end = filtered_running_periods(i, 2);
    
    fill([period_start, period_end, period_end, period_start], ...
         [0, 0, 120, 120], ...
         [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    rewards_in_period = reward_timestamps(reward_timestamps >= period_start & reward_timestamps <= period_end);

end

for i = 1:length(reward_timestamps)
    xline(reward_timestamps(i), 'r', 'LineWidth', 1.5, 'FontSize', 20);
end

yline(low_speed_threshold, '--k', 'LineWidth', 1.5, 'FontSize', 20);
yline(37.72, '--', 'LineWidth', 1.5, 'FontSize', 20, 'Color', [0.5 0.5 0.5]);


xlims = xlim;
x_offset = xlims(1) + 10;  % Same offset as above
text(x_offset, low_speed_threshold + 2, '8 cm/s', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
    'FontSize', 20, 'FontWeight', 'bold');

text(x_offset, 37.72 + 2, '37.72 cm/s', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
    'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5 0.5 0.5]);

xlabel('Time (s)', 'FontSize', 20);
ylabel('Speed (cm/s)', 'FontSize', 20);
grid on;
ylim([0, 119]);
xlim([0,480]);
set(gca, 'FontSize', 20);
hold off;

%%%%%export

% Format figure layout for better font rendering
% Prepare figure for export
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', ...
         'PaperUnits', 'Inches', ...
         'PaperSize', [pos(3), pos(4)]);

% Export to PDF without clipping
print(gcf, '2plottrace', '-dpdf', '-r0');


% Create folder if it doesn't exist
output_folder = ['running_data' RatID];
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Initialize cell array to store speed values and times
running_speeds = cell(size(filtered_running_periods, 1), 1);
running_times = cell(size(filtered_running_periods, 1), 1);

% Extract speed values within each filtered running period
for i = 1:size(filtered_running_periods, 1)
    period_start = filtered_running_periods(i, 1);
    period_end = filtered_running_periods(i, 2);

    idx = time_speed2 >= period_start & time_speed2 <= period_end;

    running_speeds{i} = speed_smth(idx);        % Speed during episode
    running_times{i} = time_speed2(idx);        % Timestamps
end

% Create final variable
TrialData = cell(1,1);
TrialData{1} = struct('times', running_times, 'speeds', running_speeds);

% Save using TrialID as the filename, variable remains 'TrialData'
save(fullfile(output_folder, [TrialID '.mat']), 'TrialData');


% Also include the full session trace
full_session.times  = time_speed2;
full_session.speeds = speed_smth;

% Append full session trace to TrialData
TrialData{2} = full_session;

% Save using TrialID as filename, variable still 'TrialData'
save(fullfile(output_folder, [TrialID '.mat']), 'TrialData');

