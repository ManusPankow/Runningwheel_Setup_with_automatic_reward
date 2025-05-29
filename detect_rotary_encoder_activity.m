

function [inactivity_start, inactivity_end] = detect_rotary_encoder_activity(pulse_times, pulse_indices, lfpSamplingRate, encoder_signal, min_pause_duration)
    % Function to detect all onsets of rotary encoder pulses and mark inactive periods
    % Outputs:
    % - inactivity_start: Start times of inactive periods
    % - inactivity_end: End times of inactive periods

    figure;
    
    % Plot rotary encoder signal
    time_vector = (1:length(encoder_signal)) / lfpSamplingRate; % Time in seconds
    plot(time_vector, encoder_signal, 'Color', [0.7 0.7 0.7]); hold on;
    
    % Mark detected pulses (onsets)
    scatter(pulse_times, ones(size(pulse_times)) * max(encoder_signal) * 0.9, ...
        50, 'b', 'filled', 'DisplayName', 'Pulse Onsets');
    
    % Detect inactive periods
    inactivity_start = [];
    inactivity_end = [];
    
    for i = 1:length(pulse_times)-1
        if (pulse_times(i+1) - pulse_times(i)) > min_pause_duration
            inactivity_start = [inactivity_start, pulse_times(i)];
            inactivity_end = [inactivity_end, pulse_times(i+1)];
        end
    end
    
    % Mark inactive sections
    for i = 1:length(inactivity_start)
        x = [inactivity_start(i), inactivity_end(i)];
        y = [min(encoder_signal), min(encoder_signal)];
        plot(x, y, 'r', 'LineWidth', 2, 'DisplayName', 'Inactive Period');
    end
    
    % Formatting
    xlabel('Time (s)');
    ylabel('Encoder Signal');
    title('Rotary Encoder Activity and Inactivity Periods');
    legend;
    grid on;
end
