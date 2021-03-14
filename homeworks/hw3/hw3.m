clear; clc; close all;
%% Question 7
variance = 1000;
M = 1000;
w = linspace(-1000*pi, 1000*pi, (2*M)+1);
m = -M:M;

noise = sqrt(variance)*randn(1,500000);

[rx, psd] = AutoCorrelation(noise,M);

figure;
stem(m, rx);
title("AutoCorrelation of White Noise");
figure;
semilogy(w,abs(psd));
title("PSD of White Noise");
ylim( [1, 3000])

%% Question 8
var = 100000;
Y = sqrt(var)*randn(10000,1);
Z = sqrt(var)*randn(10000,1);

w0 = 10000;
t = linspace(0, 10000, 10000);

X = Y.*cos(w0*t)-Z.*sin(w0*t);
[rX, psdX] = AutoCorrelation(X, 1000);

w1 = linspace(-100000*pi, 100000*pi, 2001);
figure;
semilogy(w1,abs(psdX));
title("X(t) PSD");
ylim( [1, 300000])

psdXint = (1./w1.^2).*psdX; %integration in time => 1/jw
                            %LTI => Sy(w) = |H(w)|^2 Sx(w)
figure;
semilogy(w1,abs(psdXint));
title("Integral of X(t) PSD");

%% Question 6
function [autocorr, psd] = AutoCorrelation(x, M)
    autocorr = zeros(1,(2*M)+1);
    w=linspace(-M*pi,M*pi,(2*M)+1);
    N = length(x);
    for m = 0:M
        autocorr(m+M+1) = (1/(N-m))*sum(x(1:N-m).*x(m+1:N));
    end
    for m1 = -M:-1
        autocorr(abs(m1)) = (1/(N-abs(m1)))*sum(x((abs(m1)+1):N).*x(1:N+m1));
    end
    psd = zeros(1,length(w));
    for m2 = -M:M
        psd = psd + autocorr(m2+M+1).*exp((-1j*w*m2)./(2*M+1));
    end
    
end