% Parameters
theta = pi / 3; % given angle
r_values = [0.99, 0.9, 0.7]; % different values of r
omega = linspace(-pi, pi, 1000); % frequency range

% Define the frequency response function for |H_i(e^(jω))|
Hi_magnitude = @(r, theta, omega) abs(1 - r) ./ abs(1 - 2*r*cos(theta)*exp(-1j*omega) + r^2*exp(-1j*2*omega));

% Plot magnitude response for each value of r
figure;
for i = 1:length(r_values)
    r = r_values(i);
    H_mag = Hi_magnitude(r, theta, omega); % magnitude of the frequency response
    subplot(3, 1, i);
    plot(omega, H_mag, 'LineWidth', 1.5);
    title(['Magnitude Response of |H_i(e^{j\omega})| for r = ', num2str(r)]);
    xlabel('Frequency (rad/sample)');
    ylabel('|H_i(e^{j\omega})|');
    xlim([-pi, pi]);
    grid on;
end

% Coefficients for the difference equation of Hi(z)
% Hi(z) = (1 - r) / (1 - 2*r*cos(theta)*z^(-1) + r^2*z^(-2))
for i = 1:length(r_values)
    r = r_values(i);
    b = [1 - r]; % Numerator coefficients
    a = [1, -2*r*cos(theta), r^2]; % Denominator coefficients

    % Display filter coefficients
    fprintf('Filter coefficients for r = %.2f:\n', r);
    fprintf('Numerator (b): [%.4f]\n', b);
    fprintf('Denominator (a): [%.4f, %.4f, %.4f]\n\n', a);

    % Calculate impulse response using filter coefficients
    h_n = impz(b, a, 50); % Generate 50 samples of the impulse response
    figure;
    stem(0:length(h_n)-1, h_n, 'filled');
    title(['Impulse Response h_i[n] for r = ', num2str(r)]);
    xlabel('Sample Index n');
    ylabel('h_i[n]');
    grid on;
end

% Load and play the signal
load pcm;
%sound(pcm, 8000); % Play the audio at 8 kHz

% Plot 101 samples of the signal for time indices (100:200)
figure;
plot(100:200, pcm(100:200));
title('Time-Domain Plot of pcm Signal (Indices 100:200)');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Compute and plot the magnitude of the DTFT for 1001 samples (indices 100:1100)
[X, w] = DTFT(pcm(100:1100), 1024); % Compute DTFT
figure;
plot(w, abs(X));
title('Magnitude of DTFT of pcm (Indices 100:1100)');
xlabel('Frequency (rad/sample)');
ylabel('|X(e^{j\omega})|');
xlim([-pi, pi]);
grid on;

% Find the center frequency in radians/sample (modulation frequency = 3146 Hz, fs = 8000 Hz)
theta = (2 * pi * 3146) / 8000;

% Create the IIR filter function
function y = IIRfilter(x, theta, r)
    % Initialize output y with the same length as input x
    y = zeros(size(x));
    N = length(x);
    
    % Coefficients from the transfer function Hi(z)
    b0 = 1 - r;
    a1 = -2 * r * cos(theta);
    a2 = r^2;
    
    % Apply the recursive difference equation
    for n = 1:N
        if n == 1
            y(n) = b0 * x(n); % First sample, no feedback terms
        elseif n == 2
            y(n) = b0 * x(n) - a1 * y(n-1); % Second sample, one feedback term
        else
            y(n) = b0 * x(n) - a1 * y(n-1) - a2 * y(n-2); % General case
        end
    end
end

% Apply the filter to the pcm signal
r = 0.995; % Given value of r
filtered_pcm = IIRfilter(pcm, theta, r);

% Play and plot the filtered signal for indices (100:200)
sound(filtered_pcm, 8000); % Play the filtered signal
figure;
plot(100:200, filtered_pcm(100:200));
title('Time-Domain Plot of Filtered pcm Signal (Indices 100:200)');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Compute and plot the magnitude of the DTFT for the filtered signal (indices 100:1100)
[X_filtered, w] = DTFT(filtered_pcm(100:1100), 1024);
figure;
plot(w, abs(X_filtered));
title('Magnitude of DTFT of Filtered pcm (Indices 100:1100)');
xlabel('Frequency (rad/sample)');
ylabel('|X_{filtered}(e^{j\omega})|');
xlim([-pi, pi]);
grid on;

% Zoomed plot for the range [theta - 0.02, theta + 0.02]
omega_range = (theta - 0.02):0.0001:(theta + 0.02);
X_zoomed = interp1(w, abs(X_filtered), omega_range); % Interpolate for smooth zoomed plot
figure;
plot(omega_range, X_zoomed);
title(['Zoomed Magnitude of DTFT of Filtered pcm (Around \omega = \theta)']);
xlabel('Frequency (rad/sample)');
ylabel('|X_{filtered}(e^{j\omega})|');
grid on;

% Observations:
% 1. The filter amplifies the signal near the frequency ω = θ, allowing the narrowband
%    modulated signal to be more distinct relative to background noise.
% 2. Increasing r (e.g., to 0.9999999) would create an extremely narrow peak,
%    making the filter highly selective but potentially unstable and prone to amplifying noise.
% 3. Values of r very close to 1 can create a very narrowband filter, which could
%    introduce numerical issues or instability, especially for practical applications.


% Observations:
% - As r approaches 1, the poles get closer to the unit circle, making the filter's
%   frequency response more sharply peaked at ω = ±θ. This results in a narrowband filter.
% - For lower values of r (e.g., r = 0.7), the frequency response has a broader peak,
%   indicating a wider passband.
% - The filter's nature as a bandpass filter becomes evident as it amplifies signals near
%   ω = ±θ and attenuates others.
% - The choice of r controls the selectivity of the filter: higher r values create a more
%   selective (narrower) bandpass characteristic.

% Comments:
% The filter with r values close to 1 acts as a narrowband filter, useful for extracting
% specific frequency components from a signal. As r decreases, the bandwidth of the filter
% widens, making it less selective.
