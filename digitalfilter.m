clear all;
close all;
clc;

% Filter specifications
fc = 500; % Cutoff frequency in Hz
fs = 600; % Stopband frequency in Hz
rp = 1;   % Passband ripple in dB
rs = 40;  % Stopband attenuation in dB
Fs = 2000; % Sampling frequency in Hz

% Normalize the frequencies
Wp = 2 * pi * fc / Fs; % Normalized passband frequency
Ws = 2 * pi * fs / Fs; % Normalized stopband frequency

n = 4; 

% Compute the Butterworth filter coefficients
[b, a] = butter(n, Wp / pi, 'low'); % Analog Butterworth filter

% Display the analog Butterworth filter coefficients
disp('Analog Butterworth filter coefficients:');
disp('b = '); disp(b);
disp('a = '); disp(a);

% Bilinear transformation
[bz, az] = bilinear(b, a, Fs);

% Display the digital Butterworth filter coefficients
disp('Digital Butterworth filter coefficients:');
disp('bz = '); disp(bz);
disp('az = '); disp(az);

% Frequency response of the analog filter
[H, w] = freqs(b, a, logspace(log10(Wp), log10(Ws), 1000));

% Plotting Bode plot for analog Butterworth filter
figure;
subplot(2,1,1);
semilogx(w / (2*pi), 20*log10(abs(H)));
grid on;
title('Analog Butterworth Filter Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

subplot(2,1,2);
semilogx(w / (2*pi), angle(H) * (180/pi));
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');

% Frequency response of the digital filter
[Hz, f] = freqz(bz, az, 1024, Fs);

% Plotting Bode plot for digital Butterworth filter
figure;
subplot(2,1,1);
semilogx(f, 20*log10(abs(Hz)));
grid on;
title('Digital Butterworth Filter Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

subplot(2,1,2);
semilogx(f, angle(Hz) * (180/pi));
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');

% Test signal: Sum of two sinusoids (one within passband, one within stopband)
t = 0:1/Fs:1-1/Fs;
x = sin(2*pi*100*t) + sin(2*pi*700*t); % 100 Hz and 700 Hz components

% Apply the digital filter to the test signal
y = filter(bz, az, x);

% Plot the time-domain response of the original and filtered signals together
figure;
plot(t, x, 'b'); % Original signal in blue
hold on;
plot(t, y, 'r'); % Filtered signal in red
title('Original and Filtered Combined Signals');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Signal', 'Filtered Signal');
grid on;

% Compute and plot the frequency response of the original and filtered signals
N = length(t);
X = fft(x, N);
Y = fft(y, N);
f = (0:N-1)*(Fs/N);

figure;
subplot(2,1,1);
plot(f, abs(X));
title('Frequency Response of Original Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

subplot(2,1,2);
plot(f, abs(Y));
title('Frequency Response of Filtered Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Calculate and plot the step response of the digital filter
figure;
stepz(bz, az, 1024, Fs);
title('Step Response of Digital Butterworth Filter');
grid on;

% Calculate and plot the impulse response of the digital filter
figure;
impz(bz, az, 1024, Fs);
title('Impulse Response of Digital Butterworth Filter');
grid on;

% Bode plot (magnitude and phase response) of the digital filter
figure;
[Hd, wd] = freqz(bz, az, 1024, Fs);
magnitude = 20*log10(abs(Hd));
phase = angle(Hd) * (180/pi);

subplot(2,1,1);
semilogx(wd, magnitude);
grid on;
title('Bode Plot of Digital Butterworth Filter');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

subplot(2,1,2);
semilogx(wd, phase);
grid on;
ylabel('Phase (degrees)');
xlabel('Frequency (Hz)');

% Function definition at the end of the script
function y = filter_signal(bz, az, x)
    y = filter(bz, az, x);
end
