% Define theta values
theta_values = [pi/6, pi/3, pi/2];
omega = linspace(-pi, pi, 1000); % Frequency range from -pi to pi

figure;
for i = 1:length(theta_values)
    theta = theta_values(i);
    
    % Calculate frequency response
    H_f = abs(1 - 2*cos(theta)*exp(-1j*omega) + exp(-1j*2*omega));
    
    % Plot each magnitude response
    subplot(3, 1, i);
    plot(omega, H_f);
    title(['Magnitude Response for \theta = ', num2str(theta)]);
    xlabel('\omega');
    ylabel('|H_f(e^{j\omega})|');
end

% Load and play the original signal
load nspeech1;  % Loads nspeech1 into the workspace
Fs = 8000;      % Assume a sampling rate of 8000 Hz (or adjust if specified)
%sound(nspeech1, Fs);

% Plot a segment of the original signal (samples 100 to 200)
figure;
plot(100:200, nspeech1(100:200));
title('Time Domain Plot of Original nspeech1 (Samples 100 to 200)');
xlabel('Sample Index');
ylabel('Amplitude');

% Compute DTFT of the segment (samples 100 to 1100) to find interference frequency
[X, w] = DTFT(nspeech1(100:1100), 1001); % 2048 points for better frequency resolution
figure;
plot(w, abs(X));
title('Magnitude of DTFT of Original Signal (Samples 100 to 1100)');
xlabel('Frequency (rad/sample)');
ylabel('|X(e^{j\omega})|');
xlim([-pi, pi]);

% Find the interference frequency (frequency with largest DTFT magnitude)
[~, Imax] = max(abs(X));
theta_peak = w(Imax); % Frequency of interference in radians
disp(theta_peak);

% Define FIR filter function with two zeros at e^(jθ) and e^(-jθ)
function y = FIRfilter(x, theta)
    h = [1, -2*cos(theta), 1];  % Impulse response coefficients
    y = conv(x, h, 'same');     % Apply filter via convolution
end

% Apply the filter to remove interference
filtered_signal = FIRfilter(nspeech1, theta_peak);
sound(filtered_signal, Fs);  % Listen to the filtered signal

% Plot a segment of the filtered signal (samples 100 to 200)
figure;
plot(100:200, filtered_signal(100:200));
title('Time Domain Plot of Filtered Signal (Samples 100 to 200)');
xlabel('Sample Index');
ylabel('Amplitude');

% Compute and plot DTFT of the filtered signal (samples 100 to 1100)
[X_filtered, w] = DTFT(filtered_signal(100:1100), 1001);
figure;
plot(w, abs(X_filtered));
title('Magnitude of DTFT of Filtered Signal (Samples 100 to 1100)');
xlabel('Frequency (rad/sample)');
ylabel('|H_f(e^{j\omega})|');
xlim([-pi, pi]);

% Observations:
% - The DTFT plot of the filtered signal should show reduced peaks at theta_peak,
%   indicating that the filter has attenuated the interference frequency.
