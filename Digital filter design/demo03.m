function h = LPFtrunc(N)
    omega_c = 2.0;
    n = 0:N-1;
    h = (omega_c/pi) * sinc((omega_c/pi)*(n - (N-1)/2));
end

function [X, w] = DTFT(h, Npoints)
    N = length(h);
    w = linspace(-pi, pi, Npoints);
    X = zeros(1, Npoints);
    for k = 1:Npoints
        for n = 1:N
            X(k) = X(k) + h(n) * exp(-1j * w(k) * (n-1));
        end
    end
end

% Filter sizes
N1 = 21;
N2 = 101;

% Calculate filter impulse responses
h1 = LPFtrunc(N1);
h2 = LPFtrunc(N2);

% DTFT of the filters
[X1, w1] = DTFT(h1, 512);
[X2, w2] = DTFT(h2, 512);

% Magnitude response in decibels
H1_dB = 20 * log10(abs(X1));
H2_dB = 20 * log10(abs(X2));

% Load noisy speech signal
load('nspeech2.mat');  % Assuming nspeech2 is in the .mat file

% Filter the signal using convolution
filtered1 = conv(nspeech2, h1, 'same');
filtered2 = conv(nspeech2, h2, 'same');

% Play the signals
%fs = 8000;
%sound(nspeech2, fs); % Original signal
%pause(length(nspeech2)/fs + 1);
%sound(filtered1 * 2, fs); % Filtered signal with N = 21
%pause(length(filtered1)/fs + 1);
%sound(filtered2 * 2, fs); % Filtered signal with N = 101

% Plotting the results
figure;
subplot(2,1,1);
plot(w1, abs(X1));
title('Magnitude Response of Low-Pass Filter (N = 21)');
xlabel('Frequency (\omega)');
ylabel('|H(e^{j\omega})|');
grid on;

% Marking frequency bands with color
hold on;
% Passband
hPassband = fill([-1.8, -1.8, 1.8, 1.8], [0, max(abs(X1)), max(abs(X1)), 0], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Passband

% Transition band
hLeftTransition = fill([-2.2, -2.2, -1.8, -1.8], [0, max(abs(X1)), max(abs(X1)), 0], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Left Transition band
hRightTransition = fill([1.8, 1.8, 2.2, 2.2], [0, max(abs(X1)), max(abs(X1)), 0], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Right Transition band

% Stopband
hLeftStopband = fill([-pi, -2.2, -2.2, -pi], [0, 0, max(abs(X1)), max(abs(X1))], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Left Stopband
hRightStopband = fill([2.2, 2.2, pi, pi], [0, max(abs(X1)), max(abs(X1)), 0], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Right Stopband

% Add legend
legend([hPassband, hLeftTransition, hRightTransition, hLeftStopband, hRightStopband], ...
    {'Passband', 'Left Transition Band', 'Right Transition Band', 'Left Stopband', 'Right Stopband'}, ...
    'Location', 'Best', 'AutoUpdate', 'off');

% Plot for second filter
subplot(2,1,2);
plot(w2, abs(X2));
title('Magnitude Response of Low-Pass Filter (N = 101)');
xlabel('Frequency (\omega)');
ylabel('|H(e^{j\omega})|');
grid on;

% Marking frequency bands with color
hold on;
% Passband
hPassband = fill([-1.8, -1.8, 1.8, 1.8], [0, max(abs(X1)), max(abs(X1)), 0], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Passband

% Transition band
hLeftTransition = fill([-2.2, -2.2, -1.8, -1.8], [0, max(abs(X1)), max(abs(X1)), 0], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Left Transition band
hRightTransition = fill([1.8, 1.8, 2.2, 2.2], [0, max(abs(X1)), max(abs(X1)), 0], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Right Transition band

% Stopband
hLeftStopband = fill([-pi, -2.2, -2.2, -pi], [0, 0, max(abs(X1)), max(abs(X1))], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Left Stopband
hRightStopband = fill([2.2, 2.2, pi, pi], [0, max(abs(X1)), max(abs(X1)), 0], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Right Stopband

% Add legend
legend([hPassband, hLeftTransition, hRightTransition, hLeftStopband, hRightStopband], ...
    {'Passband', 'Left Transition Band', 'Right Transition Band', 'Left Stopband', 'Right Stopband'}, ...
    'Location', 'Best', 'AutoUpdate', 'off');

% Plot in decibels
figure;
subplot(2,1,1);
plot(w1, H1_dB);
title('Magnitude Response of Low-Pass Filter (N = 21) in dB');
xlabel('Frequency (\omega)');
ylabel('Magnitude (dB)');
grid on;

subplot(2,1,2);
plot(w2, H2_dB);
title('Magnitude Response of Low-Pass Filter (N = 101) in dB');
xlabel('Frequency (\omega)');
ylabel('Magnitude (dB)');
grid on;
