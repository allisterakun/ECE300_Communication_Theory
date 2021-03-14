clc; clear; close all;

%% Setup

% initializing the signal m from .wav file and setting sampling frequence
% fs.
load handel.mat;
filename='handel.wav';
[y,Fs]=audioread(filename);
samples=[1,length(y)-(5*Fs)];
[m,fs]=audioread(filename,samples);

%sound(m, fs);
% normalizing m.
m = transpose(m(:,1) ./ max(abs(m(:,1))));

N = length(m);
T = N/fs;
fc = 55000;

%upsampling
scale = 3*ceil(fc/fs);
N_up = scale * N;
fs_up = scale * fs;
t = linspace(0,T,N);
t_up = linspace(0,T,N_up);
wd = linspace(-pi, pi, N_up);
f = wd *fs_up / (2 * pi);
fdemod = linspace(-pi, pi, N) * fs / (2 * pi);
fc = 55000;

sSin = sin(2*pi * fc * t_up);
m_P = rms(m).^2;
W_m = obw(m, fs);
%% FM Mod & Demod
Ac_FM = 1.21;
k_FM = 70000;
m_FM = FM_mod(m, fs, Ac_FM, 55000, k_FM);
m_FMfreq = fftshift(fft(m_FM)/fs_up);

figure('Position', [0 0 1024 1024]);
subplot(2,1,1)
plot(t_up,m_FM);
title("FM Modulated m(t)");xlabel("Time (s)"); ylabel("u(t)");
subplot(2,1,2)
semilogy(f, abs(m_FMfreq));
title("Fourier Transform of FM Modulated m(t)");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
suptitle("FM Modulation");

m_FM_demod = FM_demod(m_FM, fs_up, 2000,f, scale);
m_FM_demodfreq = fftshift(fft(m_FM_demod)/fs);

figure('Position', [0 0 1024 1024]);
subplot(2,1,1)
plot(t,m_FM_demod);
title("FM Demodulated m(t)");xlabel("Time (s)"); ylabel("u(t)");
subplot(2,1,2)
semilogy(fdemod, abs(m_FM_demodfreq));
title("Fourier Transform of FM Demodulated m(t)");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
suptitle("FM Demodulation");

%sound(m_FM_demod, fs);

m_FM_P = rms(m_FM_demod).^2;
W_FM = obw(m_FM_demod, fs);

var = [0.01, 0.05, 0.1];
FM_noiseless = rms(m_FM_demod).^2;

for i = 1:3

    v = var(i);
    noise = sqrt(var(i)) * randn(1,length(m_FM));
    m_FM_noise = m_FM + noise;
    m_FM_noisefreq = fftshift(fft(m_FM_noise)/fs_up);
    m_FM_noise_demod = FM_demod(m_FM_noise, fs_up, 2000,f, scale);
    m_FM_noise_demod = lowpass(m_FM_noise_demod, 2000, fs);
    m_FM_noise_demodfreq = fftshift(fft(m_FM_noise_demod)/fs);
    
    figure('Position', [0 0 1024 1024]);
    subplot(2, 2, 1);
    semilogy(f, abs(m_FMfreq));
    title("  Magnitude of Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    subplot(2, 2, 2);
    semilogy(fdemod, abs(m_FM_demodfreq));
    title("  Magnitude of Modulated Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    subplot(2, 2, 3);
    semilogy(f, abs(m_FM_noisefreq));
    title("  Magnitude of Modulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    subplot(2, 2, 4);
    semilogy(fdemod, abs(m_FM_noise_demodfreq));
    title("  Magnitude of Demodulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    suptitle("FM with Noise \sigma^2 = " + var(i));
    
    noisy = rms(m_FM_noise_demod - m_FM_demod).^2;
    snr_FM = 10 * log10(FM_noiseless/noisy);
    display("SNR for FM with variance" + var(i) + ": " + snr_FM);
    
    snr_FM_theo = (Ac_FM^2 * k_FM^2 * m_FM_P) / (4 * var(i) * W_FM^3 / fs);
    display("Theoretical SNR for FM with variance " + var(i) + ": " + snr_FM_theo);
end

figure('Position', [0 0 1024 1024]);
for j = 1:3
    var = 0.1;
    Ac_FM = 1.21;
    k_FM = 70000 * 2^(j-2);
    m_FM = FM_mod(m, fs, Ac_FM, 55000, k_FM);
    m_FMfreq = fftshift(fft(m_FM)/fs_up);
    m_FM_demod = FM_demod(m_FM, fs_up, 2000,f, scale);
    m_FM_demodfreq = fftshift(fft(m_FM_demod)/fs);
    FM_noiseless = rms(m_FM_demod).^2;
    noise = sqrt(var) * randn(1,length(m_FM));
    m_FM_noise = m_FM + noise;
    m_FM_noisefreq = fftshift(fft(m_FM_noise)/fs);
    m_FM_noise_demod = PM_demod(m_FM_noise, 1, fs_up, 2000, scale, sSin);
    m_FM_noise_demod = lowpass(m_FM_noise_demod, 2000, fs);
    m_FM_noise_demodfreq = fftshift(fft(m_FM_noise_demod)/fs);
    
    noisy = rms(m_FM_noise_demod - m_FM_demod).^2;
    snr_FM = 10 * log10(FM_noiseless/noisy);
    display("SNR for FM with k = " + k_FM + "(variance = "+var+") is: " + snr_FM);
    
    snr_FM_theo = (Ac_FM^2 * k_FM^2 * m_FM_P) / (4 * var * W_FM^3 / fs);
    display("Theoretical SNR for FM with k = " + k_FM + "(variance = "+var+") is: " + snr_FM_theo);
    
    subplot(2, 2, 1);
    str = strcat('k = ' , int2str(k_FM));
    semilogy(f, abs(m_FMfreq), 'DisplayName', str);
    title("  Magnitude of Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    legend;
    hold on
    
    subplot(2, 2, 2);
    semilogy(fdemod, abs(m_FM_demodfreq), 'DisplayName', str);
    title("  Magnitude of Modulated Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    legend;
    hold on
    
    subplot(2, 2, 3);
    semilogy(f, abs(m_FM_noisefreq), 'DisplayName', str);
    title("  Magnitude of Modulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    legend;
    hold on
    
    subplot(2, 2, 4);
    semilogy(fdemod, abs(m_FM_noise_demodfreq), 'DisplayName', str);
    title("  Magnitude of Demodulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    legend;
    hold on;
    
end
suptitle("FM with Noise \sigma^2 = " + var);
hold off;

%% PM Mod & Demod
Ac_PM = 0.6019;
k_PM = 2;
m_PM = PM_mod(m, fs, Ac_PM, 55000, k_PM);
m_PMfreq = fftshift(fft(m_PM)/fs_up);

figure('Position', [0 0 1024 1024]);
subplot(2,1,1)
plot(t_up,m_PM);
title("PM Modulated m(t)");xlabel("Time (s)"); ylabel("u(t)");
subplot(2,1,2)
semilogy(f, abs(m_PMfreq));
title("Fourier Transform of PM Modulated m(t)");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
suptitle("PM Modulation");

m_PM_demod = PM_demod(m_PM, 2, fs_up, 2000, scale, sSin);
m_PM_demodfreq = fftshift(fft(m_PM_demod)/fs);

figure('Position', [0 0 1024 1024]);
subplot(2,1,1)
plot(t,m_PM_demod);
title("PM Demodulated m(t)");xlabel("Time (s)"); ylabel("u(t)");
subplot(2,1,2)
semilogy(fdemod, abs(m_PM_demodfreq));
title("Fourier Transform of PM Demodulated m(t)");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
suptitle("PM Demodulation");

%sound(m_PM_demod, fs);
m_PM_P = rms(m_PM_demod).^2;
W_PM = obw(m_PM_demod, fs);

var = [0.01, 0.05, 0.1];
PM_noiseless = rms(m_PM_demod).^2;

for i = 1:3

    v = var(i);
    noise = sqrt(var(i)) * randn(1,length(m_PM));
    m_PM_noise = m_PM + noise;
    m_PM_noisefreq = fftshift(fft(m_PM_noise)/fs);
    m_PM_noise_demod = PM_demod(m_PM_noise, 1, fs_up, 2000, scale, sSin);
    m_PM_noise_demod = lowpass(m_PM_noise_demod, 2000, fs);
    m_PM_noise_demodfreq = fftshift(fft(m_PM_noise_demod)/fs);
    
    figure('Position', [0 0 1024 1024]);
    subplot(2, 2, 1);
    semilogy(f, abs(m_PMfreq));
    title("  Magnitude of Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    subplot(2, 2, 2);
    semilogy(fdemod, abs(m_PM_demodfreq));
    title("  Magnitude of Modulated Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    subplot(2, 2, 3);
    semilogy(f, abs(m_PM_noisefreq));
    title("  Magnitude of Modulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    subplot(2, 2, 4);
    semilogy(fdemod, abs(m_PM_noise_demodfreq));
    title("  Magnitude of Demodulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    suptitle("PM with Noise \sigma^2 = " + var(i));
    
    noisy = rms(m_PM_noise_demod - m_PM_demod).^2;
    snr_PM = 10 * log10(PM_noiseless/noisy);
    display("SNR for PM with variance" + var(i) + ": " + snr_PM);
    
    snr_PM_theo = (Ac_PM^2 * k_PM^2 * m_PM_P) / (4 * var(i) * W_PM / fs);
    display("Theoretical SNR for PM with variance " + var(i) + ": " + snr_PM_theo);
end

figure('Position', [0 0 1024 1024]);
for j = 1:3
    var = 0.1;
    Ac_PM = 10827/10000;
    k_PM = 2^(j-1);
    m_PM = PM_mod(m, fs, Ac_PM, 55000, k_PM);
    m_PMfreq = fftshift(fft(m_PM)/fs_up);
    m_PM_demod = PM_demod(m_PM, 2, fs_up, 2000, scale, sSin);
    m_PM_demodfreq = fftshift(fft(m_PM_demod)/fs);
    PM_noiseless = rms(m_PM_demod).^2;
    noise = sqrt(var) * randn(1,length(m_PM));
    m_PM_noise = m_PM + noise;
    m_PM_noisefreq = fftshift(fft(m_PM_noise)/fs);
    m_PM_noise_demod = PM_demod(m_PM_noise, 1, fs_up, 2000, scale, sSin);
    m_PM_noise_demod = lowpass(m_PM_noise_demod, 2000, fs);
    m_PM_noise_demodfreq = fftshift(fft(m_PM_noise_demod)/fs);
    
    noisy = rms(m_PM_noise_demod - m_PM_demod).^2;
    snr_PM = 10 * log10(PM_noiseless/noisy);
    display("SNR for PM with k = " + k_PM + "(variance = "+var+") is: " + snr_PM);
    
    snr_PM_theo = (Ac_PM^2 * k_PM^2 * m_PM_P) / (4 * var * W_PM / fs);
    display("Theoretical SNR for PM with k = " + k_PM + "(variance = "+var+") is: " + snr_PM_theo);
    
    subplot(2, 2, 1);
    str = strcat('k = ' , int2str(k_PM));
    semilogy(f, abs(m_PMfreq), 'DisplayName', str);
    title("  Magnitude of Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    legend;
    hold on
    
    subplot(2, 2, 2);
    semilogy(fdemod, abs(m_PM_demodfreq), 'DisplayName', str);
    title("  Magnitude of Modulated Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    legend;
    hold on
    
    subplot(2, 2, 3);
    semilogy(f, abs(m_PM_noisefreq), 'DisplayName', str);
    title("  Magnitude of Modulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    legend;
    hold on
    
    subplot(2, 2, 4);
    semilogy(fdemod, abs(m_PM_noise_demodfreq), 'DisplayName', str);
    title("  Magnitude of Demodulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    legend;
    hold on;
    
end
suptitle("PM with Noise \sigma^2 = " + var);
hold off;
%% AM Mod & Demod
Ac_AM = 87/250;
a = 2;

m_AM = AM_mod(m, fs, Ac_AM, 55000,a);
m_AMfreq = fftshift(fft(m_AM)/fs_up);

figure('Position', [0 0 1024 1024]);
subplot(2,1,1)
plot(t_up,m_AM);
title("AM Modulated m(t)");xlabel("Time (s)"); ylabel("u(t)");
subplot(2,1,2)
semilogy(f, abs(m_AMfreq));
title("Fourier Transform of AM Modulated m(t)");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
suptitle("AM Modulation");

m_AM_demod = AM_demod(m_AM, fs_up, 2000, scale);
m_AM_demodfreq = fftshift(fft(m_AM_demod)/fs);

figure('Position', [0 0 1024 1024]);
subplot(2,1,1)
plot(t,m_AM_demod);
title("AM Demodulated m(t)");xlabel("Time (s)"); ylabel("u(t)");
subplot(2,1,2)
semilogy(fdemod, abs(m_AM_demodfreq));
title("Fourier Transform of AM Demodulated m(t)");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
suptitle("AM Demodulation");

%sound(m_AM_demod, fs);
m_AM_P = rms(m_AM_demod).^2;
W_AM = obw(m_AM_demod, fs);

var = [0.01, 0.05, 0.1];
AM_noiseless = rms(m_AM_demod).^2;

for i = 1:3

    v = var(i);
    noise = sqrt(var(i)) * randn(1,length(m_AM));
    m_AM_noise = m_AM + noise;
    m_AM_noisefreq = fftshift(fft(m_AM_noise)/fs);
    m_AM_noise_demod = AM_demod(m_AM_noise, fs_up, 2000, scale);
    m_AM_noise_demod = lowpass(m_AM_noise_demod, 2000, fs);
    m_AM_noise_demodfreq = fftshift(fft(m_AM_noise_demod)/fs);
    figure('Position', [0 0 1024 1024]);
    subplot(2, 2, 1);
    semilogy(f, abs(m_AMfreq));
    title("  Magnitude of Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    subplot(2, 2, 2);
    semilogy(fdemod, abs(m_AM_demodfreq));
    title("  Magnitude of Modulated Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    subplot(2, 2, 3);
    semilogy(f, abs(m_AM_noisefreq));
    title("  Magnitude of Modulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    subplot(2, 2, 4);
    semilogy(fdemod, abs(m_AM_noise_demodfreq));
    title("  Magnitude of Demodulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    suptitle("AM with Noise \sigma^2 = " + var(i));
    
    noisy = rms(m_AM_noise_demod - m_AM_demod).^2;
    snr_AM = 10 * log10(AM_noiseless/noisy);
    display("SNR for AM with variance" + var(i) + ": " + snr_AM);
    
    snr_AM_theo = (Ac_AM^2 * a^2 * m_AM_P) / (4 * var(i) * W_AM / fs);
    display("Theoretical SNR for AM with variance " + var(i) + ": " + snr_AM_theo);
end
%sound(m_AM_noise_demod, fs);
figure('Position', [0 0 1024 1024]);
for j = 1:3
    var = 0.1;
    Ac_AM = 87/250;
    a = 2^(j-1);

    m_AM = AM_mod(m, fs, Ac_AM, 55000,a);
    m_AMfreq = fftshift(fft(m_AM)/fs_up);
    m_AM_demod = AM_demod(m_AM, fs_up, 2000, scale);
    m_AM_demodfreq = fftshift(fft(m_AM_demod)/fs);
    AM_noiseless = rms(m_AM_demod).^2;
    noise = sqrt(var) * randn(1,length(m_AM));
    m_AM_noise = m_AM + noise;
    m_AM_noisefreq = fftshift(fft(m_AM_noise)/fs);
    m_AM_noise_demod = AM_demod(m_AM_noise, fs_up, 2000, scale);
    m_AM_noise_demod = lowpass(m_AM_noise_demod, 2000, fs);
    m_AM_noise_demodfreq = fftshift(fft(m_AM_noise_demod)/fs);
    
    noisy = rms(m_AM_noise_demod - m_AM_demod).^2;
    snr_AM = 10 * log10(AM_noiseless/noisy);
    display("SNR for AM with a = " + a + "(variance = "+var+") is: " + snr_AM);
    
    snr_AM_theo = (Ac_AM^2 * a^2 * m_AM_P) / (4 * var * W_AM / fs);
    display("Theoretical SNR for AM with a = " + a + "(variance = "+var+") is: " + snr_AM_theo);
    
    subplot(2, 2, 1);
    str = strcat('a = ' , int2str(a));
    semilogy(f, abs(m_AMfreq), 'DisplayName', str);
    title("  Magnitude of Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    legend;
    hold on
    
    subplot(2, 2, 2);
    semilogy(fdemod, abs(m_AM_demodfreq), 'DisplayName', str);
    title("  Magnitude of Modulated Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    legend;
    hold on
    
    subplot(2, 2, 3);
    semilogy(f, abs(m_AM_noisefreq), 'DisplayName', str);
    title("  Magnitude of Modulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    legend;
    hold on
    
    subplot(2, 2, 4);
    semilogy(fdemod, abs(m_AM_noise_demodfreq), 'DisplayName', str);
    title("  Magnitude of Demodulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    legend;
    hold on;
    
end
suptitle("AM with Noise \sigma^2 = " + var);
hold off;
%% SSB Mod & Demod
Ac_SSB = 262/125;

m_SSB = SSB_mod(m, fs, Ac_SSB, 55000);
m_SSBfreq = fftshift(fft(m_SSB)/fs_up);

figure('Position', [0 0 1024 1024]);
subplot(2,1,1)
plot(t_up,m_SSB);
title("SSB Modulated m(t)");xlabel("Time (s)"); ylabel("u(t)");
subplot(2,1,2)
semilogy(f, abs(m_SSBfreq));
title("Fourier Transform of SSB Modulated m(t)");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
suptitle("SSB Modulation");

m_SSB_demod = SSB_demod(m_SSB, fs_up, 2000, scale, sSin);
m_SSB_demodfreq = fftshift(fft(m_SSB_demod)/fs);

figure('Position', [0 0 1024 1024]);
subplot(2,1,1)
plot(t,m_SSB_demod);
title("SSB Demodulated m(t)");xlabel("Time (s)"); ylabel("u(t)");
subplot(2,1,2)
semilogy(fdemod, abs(m_SSB_demodfreq));
title("Fourier Transform of SSB Demodulated m(t)");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
suptitle("SSB Demodulation");

%sound(m_SSB_demod, fs);

m_SSB_P = rms(m_SSB_demod).^2;
W_SSB = obw(m_SSB_demod, fs);

var = [0.01, 0.05, 0.1];
SSB_noiseless = rms(m_SSB_demod).^2;

for i = 1:3

    v = var(i);
    noise = sqrt(var(i)) * randn(1,length(m_SSB));
    m_SSB_noise = m_SSB + noise;
    m_SSB_noisefreq = fftshift(fft(m_SSB_noise)/fs);
    m_SSB_noise_demod = SSB_demod(m_SSB_noise, fs_up, 2000, scale, sSin);
    m_SSB_noise_demod = lowpass(m_SSB_noise_demod, 2000, fs);
    m_SSB_noise_demodfreq = fftshift(fft(m_SSB_noise_demod)/fs);
    
    figure('Position', [0 0 1024 1024]);
    subplot(2, 2, 1);
    semilogy(f, abs(m_SSBfreq));
    title("  Magnitude of Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    subplot(2, 2, 2);
    semilogy(fdemod, abs(m_SSB_demodfreq));
    title("  Magnitude of Modulated Noise Free Input");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    subplot(2, 2, 3);
    semilogy(f, abs(m_SSB_noisefreq));
    title("  Magnitude of Modulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    subplot(2, 2, 4);
    semilogy(fdemod, abs(m_SSB_noise_demodfreq));
    title("  Magnitude of Demodulated Noisy Signal");
    xlabel("Frequency (Hz)"); ylabel("Amplitude");
    suptitle("SSB with Noise \sigma^2 = " + var(i));
    
    noisy = rms(m_SSB_noise_demod - m_SSB_demod).^2;
    snr_SSB = 10 * log10(SSB_noiseless/noisy);
    display("SNR for SSB with variance" + var(i) + ": " + snr_SSB);
    
    snr_SSB_theo = (Ac_SSB^2 * m_SSB_P) / (2 * var(i) * W_SSB / fs);
    display("Theoretical SNR for SSB with variance " + var(i) + ": " + snr_SSB_theo);
end

%% Functions
%%%%%%%%%%%%%%%%%%%% FM
function [output] = FM_mod(m, fs, Ac, fc, k)
    N = length(m);
    T = N/fs;
    t = linspace(0,T,N);
    wc = fc * 2 * pi;
    
    %upsampling
    scale = 3*ceil(fc/fs);
    N_up = scale * N;
    fs_up = scale * fs;
    t_up = linspace(0,T,N_up);
    m_up = interp1(t,m,t_up);
    
    
    theta = cumsum(m_up)/fs_up * 2 * pi * k;
    output = (Ac * cos(wc * t_up + theta));
    
end

function [output] = FM_demod(m, fs, fcut, f, scale)
    %m = abs(m);
    w = 2 * pi * f;
    output = real(ifft(((fft(m)/fs)) .* (1j * w)));
    output(output<0) = 0;
    output = lowpass(output,fcut,fs);
    output = downsample(output,scale);
end

%%%%%%%%%%%%%%%%%%%% PM
function [output] = PM_mod(m, fs, Ac, fc, k)
    N = length(m);
    T = N/fs;
    t = linspace(0,T,N);
    wc = fc * 2 * pi;
    
    %upsampling
    scale = 3*ceil(fc/fs);
    N_up = scale * N;
    %fs_up = scale * fs;
    t_up = linspace(0,T,N_up);
    m_up = interp1(t,m,t_up);
    
    theta = m_up * k;
    output = (Ac * cos(wc * t_up + theta));
end

function [output] = PM_demod(m, Ac, fs, fcut, scale, sSin)
    output = Ac .* m .* sSin;
    output = lowpass(output, fcut, fs);
    output = downsample(output, scale);
    
end

%%%%%%%%%%%%%%%%%%%% AM
function [output] = AM_mod(m, fs, Ac, fc,a)
    N = length(m);
    T = N/fs;
    t = linspace(0,T,N);
    wc = fc * 2 * pi;
    
    %upsampling
    scale = 3*ceil(fc/fs);
    N_up = scale * N;
    %fs_up = scale * fs;
    t_up = linspace(0,T,N_up);
    m_up = interp1(t,m,t_up);
    
    cSig = cos(wc * t_up);
    output = Ac * (1+ a*m_up) .* cSig;
    %output = highpass(output, 1000, fs);
end

function [output] = AM_demod(m, fs, fcut, scale)
    
    %m = m+0.5;
    %m(m<0) = 0;
    m = abs(m);
    %t_up = linspace(0, length(m)/(fs/scale), length(m)*scale);
    %output = cos(2*pi*55000.*t_up).*m;
    output = lowpass(m,fcut,fs);
    %output = highpass(output, 2, fs);
    output = downsample(output, scale);
    %output = lowpass(output, 20000, fs);
end

%%%%%%%%%%%%%%%%%%%% SSB
function [output] = SSB_mod(m, fs, Ac, fc)
    N = length(m);
    T = N/fs;
    t = linspace(0,T,N);
    wc = fc * 2 * pi;
    
    %upsampling
    scale = 3*ceil(fc/fs);
    N_up = scale * N;
    fs_up = scale * fs;
    t_up = linspace(0,T,N_up);
    m_up = interp1(t,m,t_up);
    wd = linspace(-pi, pi, N_up);
    f = wd *fs_up / (2 * pi);
    
    cSig = cos(wc * t_up);
    sSig = sin(wc * t_up);
    mH = real(ifft(fft(m_up).*(-1j.*sign(f))));
    output = Ac * (m_up .* cSig - mH .* sSig);
    %output = highpass(output, 5000, fs);
end

function [output] = SSB_demod(m, fs, fcut, scale, sSin)
    m = m .* sSin;
    output = lowpass(m,fcut,fs);
    output = downsample(output, scale);
end