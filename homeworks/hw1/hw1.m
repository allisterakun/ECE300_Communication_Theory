close all; clear; clc;

%% Quesiton 5
N = 10000;
t = linspace(-10, 10, N); fs = N/20;
frequency = 1; T = 1/frequency;
rect = t;
rect (rect<0|rect>T) = 0; rect (rect>0&rect<T) = 1;
x = rect .* cos(2*pi*frequency*t);
X = (fft(x.*rect))/fs;
Xs = fftshift(X);
normofX = (sum(abs(Xs.^2)))*(500/N);
normofx = (sum(abs(x.^2)))*(20/N);

%% Question 6
wd = linspace(-pi, pi, N); f = wd*fs/(2*pi);
y = abs(circshift(x, -3*fs));
X = fft(x); Xs = fftshift(X);
Y = fft(y); Ys = fftshift(Y);

figure;
subplot(2,2,1);
plot(f, abs(Xs));
title("X(f)");
xlabel("f (Hz)");
xlim([-50 50])


subplot(2,2,3);
plot(f, abs(Ys));
title("Y(f)");
xlabel("f (Hz)");
xlim([-50 50])


subplot(2,2,2);
plot(f, angle(Xs));
title("Phase of X(f)");
xlabel("f (Hz)");
xlim([-50 50])


subplot(2,2,4);
plot(f, angle(Ys));
title("Phase of Y(f)");
xlabel("f (Hz)");
xlim([-50 50])


%% Question 7

theta = pi/3;
zbp = y.*cos((64*pi*t)+theta); 
zbb = y.*(cos(theta)+1j*sin(theta));

Z = fft(zbp);
Zs = fftshift(Z);

figure;
subplot(2,3,1);
plot(t, zbp);
title("Bandpass Signal zbp(t)");
xlabel("t");
xlim([-3.5 -1.5]);

subplot(2,3,4);
plot(t, abs(zbb));
title("Baseband Signal zbb(t)");
xlabel("t");
xlim([-3.5 -1.5]);

subplot(2,3,2);
plot(f, abs(Ys));
title("Y(w)");
xlabel("Frequency (Hz)");
xlim([-50, 50]);

subplot(2,3,5);
plot(f, abs(Zs));
title("Z(w)");
xlabel("Frequency (Hz)");
xlim([-50, 50]);

subplot(2,3,3);
plot(f, angle(Ys));
title("Phase of Y(w)");
xlabel("Frequency (Hz)");


subplot(2,3,6);
plot(f, angle(Zs));
title("Phase of Z(w)");
xlabel("Frequency (Hz)");


%% Question 8
[Hx, HX] = HilbertTransform(x, f);
[Hy, HY] = HilbertTransform(y, f);
[Hz, HZ] = HilbertTransform(zbp, f);

figure;
subplot(2,3,1);
plot(f, abs(Hx));
title("H(w) of x");
xlabel("Frequency (Hz)")

subplot(2,3,2);
plot(f, abs(Hy));
title("H(w) of y");
xlabel("Frequency (Hz)")

subplot(2,3,3);
plot(f, abs(Hz));
title("H(w) of z");
xlabel("Frequency (Hz)")

subplot(2,3,4);
plot(t, HX);
title("H(t) of x");
xlabel("t")
xlim([-10, 10]);

subplot(2,3,5);
plot(t, HY);
title("H(t) of y");
xlabel("t")
xlim([-10, 10]);

subplot(2,3,6);
plot(t, HZ);
title("H(t) of z");
xlabel("t")
xlim([-10, 10]);

%% Question 9
orthox = dot(x, HX);
orthoy = dot(y, HY);
orthoz = dot(zbp, HZ);

% The dot products of x, y, z, with their Hilbert transforms are very close to 0
% which tells us that they are orthogonal.
%% Functions
function [outf, outt] = HilbertTransform(signal, f)
    outf1 = sign(f) .* -1i .* fftshift(fft(signal));
    outf = real(outf1);
    outt = real(ifft(fftshift(outf1)));
end
