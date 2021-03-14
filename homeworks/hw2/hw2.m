close all; clear; clc;

%% Setup

[m, fs] = audioread("car_x.wav");
mt= transpose(m(:,1) ./ max(abs(m(:,1))));          %normalizing m.
N = length(m); T = N / fs;
scale = 20;
% fs is approximately 11k, so 2*(550k/11k) = 10,
% give it some room, we make it 20.
fsScaled = fs * scale; Nscaled = N * scale;         %upsacling N and fs.
t = linspace(0, T, N); 
tScaled = linspace(0, T, Nscaled);
mtScaled = interp1(t, mt, tScaled);

M = fftshift(fft(mtScaled)/fsScaled);               %Fourier trans
wd = linspace(-pi, pi, Nscaled);
f = wd * fsScaled / (2*pi);

%% Question 1
figure;

subplot(3, 2, [1, 2]);
plot(tScaled, mtScaled);
title("m(t) in time");
xlabel("time (s)"); 

subplot(3, 2, 3);
semilogy(f, abs(M));
title("Amplitude of M(w) - Semilog");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
xlim([-25000, 25000]);

subplot(3, 2, 4);
plot(f, abs(M));
title("Amplitude of M(w)");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
xlim([-25000, 25000]);

subplot(3, 2, [5, 6]);
plot(f, unwrap(angle(M)));
title("Angle of M(w)"); ylabel("Phase");
xlabel("Frequency (Hz)");

%By looking at the graph of M(w), the bandwidth is approximately 6000 Hz.

%% Question 2

fc = 55000;
carrierSignal = cos(2*pi * fc * tScaled);

%% Question 3

m_dsbsc = mtScaled .* carrierSignal;
m_dsb = mtScaled .* carrierSignal + carrierSignal;

M_dsbsc = fftshift(fft(m_dsbsc) / fsScaled);
M_dsb = fftshift(fft(m_dsb) / fsScaled);

figure;

subplot(2, 4, [1, 2]);
plot(tScaled, m_dsbsc);
title("DSB-SC in time domain");
xlabel("Time (s)");

subplot(2, 4, [3, 4]);
plot(tScaled, m_dsb);
title("DSB in time domain");
xlabel("Time (s)");

subplot(2, 4, 5);
semilogy(f, abs(M_dsbsc));
title("DSB-SC in frequency domain");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
%xlim([-80000, 80000]);

subplot(2, 4, 7);
semilogy(f, abs(M_dsb));
title("DSB in frequency domain");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
%xlim([-80000, 80000]); ylim([0 0.01]);

subplot(2, 4, 6)
plot(f, unwrap(angle(M_dsbsc)));
title("Angle of DSB-SC"); ylabel("Phase");
xlabel("Frequency (Hz)");

subplot(2, 4, 8)
plot(f, unwrap(angle(M_dsb)));
title("Angle of DSB"); ylabel("Phase");
xlabel("Frequency (Hz)");

%Comparing to DSB-SC, the DSB signal has an additional cosine 
%added to it.

%% Question 4

sins = sin(2*pi * fc * tScaled);

mH = real(ifft(fft(mtScaled).*(-1j.*sign(f))));

lowerssb = mtScaled .*carrierSignal + mH .* sins;
upperssb = mtScaled .*carrierSignal - mH .* sins;

LowerSSB = fftshift(fft(lowerssb)/fsScaled);
UpperSSB = fftshift(fft(upperssb)/fsScaled);


figure;
subplot(2,4,[1, 2]);
plot(tScaled, lowerssb);
title("Lower SSB in time domain");
xlabel("Time (s)");

subplot(2,4,[3 ,4]);
plot(tScaled, upperssb);
title("Upper SSB in time domain");
xlabel("Time (s)");

subplot(2,4,[5, 6]);
semilogy(f, abs(LowerSSB));
title("Lower SSB in frequency domain");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
%xlim([-100000, 100000]);

subplot(2,4,[7, 8]);
semilogy(f, abs(UpperSSB));
title("Upper SSB in frequency domain");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
%xlim([-100000, 100000]);

%% Question 5

A = 1;
conventional = (1+mtScaled) .* carrierSignal * A;
CONVENTIONAL = fftshift(fft(conventional)/fsScaled);

figure;
subplot(2,1,1);
plot(tScaled, conventional);
title("Conventional Signal in time domain");
xlabel("Time (s)");

subplot(2,1,2);
semilogy(f, abs(CONVENTIONAL));
title("Conventional Signal in frequency domain");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
%xlim([-80000, 80000]); ylim([0 0.4]);

%% Question 6

conventionalRect = abs(conventional);
CONVENTIONALrect = fftshift(fft(conventionalRect)/fsScaled);

figure;
subplot(2,1,1);
plot(tScaled, conventionalRect);
title("Rectified Conventional Signal in time domain");
xlabel("Time (s)");

subplot(2,1,2);
semilogy(f, abs(CONVENTIONALrect));
title("Rectified Conventional Signal in frequency domain");
xlabel("Frequency (Hz)"); ylabel("Amplitude");
%xlim([-50000, 50000]); ylim([0 0.01]);

conventionalRectLow = lowpass(conventionalRect, 2500, fsScaled);
%conventionalRectLow = lowpass(conventionalRectLow, 2500, fsScaled);

figure;
subplot(2,1,1);
plot(tScaled, conventionalRectLow);
title("Low-passed Signal in time domain");
xlabel("Time (s)");

subplot(2,1,2);
semilogy(f, abs(CONVENTIONALrect));
title("Low-passed Signal in frequency domain");
xlabel("Frequency (Hz)");ylabel("Amplitude");
%xlim([-50000, 50000]); ylim([0 0.01]);

%play the sound
conventionalRectPlay = downsample(conventionalRectLow, scale);
%figure;
%plot(t, conventionalRectPlay);
sound(conventionalRectPlay, fs);

