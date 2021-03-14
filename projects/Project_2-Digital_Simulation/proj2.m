close all; clc; clear;

nExample = 100000; % Number of symbols

%% Question 1
constellation = [ 1, -1, 1j, -1j ];
var = 0.1;

% generating symbol
symbol = constellation(randsample(length(constellation),nExample,'true'));

% generating AWGN
noise = sqrt(var)*(randn(1,nExample) + 1j*randn(1,nExample));

% adding noise to symbol to get noisy symbol
symbolNoisy = symbol + noise;

% estimated symbol vector based on least squares decision method
LSestimates = LSestimation(constellation, symbolNoisy);

%% Question 2,3,4
N0 = 5;
SNRpBit = 15;
[nf, received, sf] = scale(constellation, symbol, N0, SNRpBit, false);
err = calcError(nf, received);
[nfDiff, receivedDiff, sfDiff] = scale(constellation, symbol, N0, SNRpBit, true);
errDiff = calcError(nfDiff, receivedDiff);

%% Question 5
N0 = 1;
SNRvalue = -4:20;

% BINARY ANTIPODAL
    
    tic;
    % constelltion for binary antipodal
    constellationBinAp = [1j, -1j];
    
    % generate symbols for binary antipodl
    symbolBinAp = symbolGen(constellationBinAp, nExample);
    
    % initializing theoretical and actual binary antipodal Perror vector
    PerrorTheoBinAp = zeros(1, length(SNRvalue));
    PerrorActuBinAp = zeros(1, length(SNRvalue));
    
    M = length(constellationBinAp);
    
    for i = 1:length(SNRvalue)
        % simulation
        [nfBinAp,receivedBinAp, sfAp] = scale(constellationBinAp,symbolBinAp,N0,i,false);
        
        % actual error
        PerrorActuBinAp(i) = calcError(nfBinAp, receivedBinAp);
        
        % theoretical error - page 406 from the book
        PerrorTheoBinAp(i) = qfunc(sfAp * sqrt(4/N0));
    end
    
    fprintf('Time taken for binary antipodal with N = %d: %.2f s\n', nExample, toc());
    
% BINARY ORTHOGONAL
    tic;
    % constelltion for binary orthogonal
    constellationBinOr = [1, 1j];
    
    % generate symbols for binary antipodl
    symbolBinOr = symbolGen(constellationBinOr, nExample);
    
    % initializing theoretical and actual binary antipodal Perror vector
    PerrorTheoBinOr = zeros(1, length(SNRvalue));
    PerrorActuBinOr = zeros(1, length(SNRvalue));
    
    M = length(constellationBinOr);
    
    for i = 1:length(SNRvalue)
        % energy
        Eng = 10^(i/20) * N0/2 * log2(M);
        
        % simulation
        [nfBinOr,receivedBinOr,sfOr] = scale(constellationBinOr,symbolBinOr,N0,i,false);
        
        % actual error
        PerrorActuBinOr(i) = calcError(nfBinOr, receivedBinOr);
        
        % theoretical error - page 408 from the book
        PerrorTheoBinOr(i) = qfunc(sfOr * sqrt(2/N0));
    end
    
    fprintf('Time taken for binary orthogonal with N = %d: %.2f s\n', nExample, toc());
% PLOTTING GRAPH

    %fprintf('plotting\n');
    figure();
    hold on;
    scatter(SNRvalue,PerrorActuBinAp,'.','b','DisplayName','Antipodal');
    plot(SNRvalue,PerrorTheoBinAp,'Color','blue','DisplayName','Theoretical Antipodal');
    scatter(SNRvalue,PerrorActuBinOr,'.','r','DisplayName','Orthogonal');
    plot(SNRvalue,PerrorTheoBinOr,'Color','red','DisplayName','Theoretical Orthogonal');
    hold off;
    title('Binary Perror vs. SNR');
    ylabel('P_{error}');
    xlabel('SNR/bit (dB)');
    legend();
    
%% Question 6
% PSK
    tic;
    SNRvalue = -4:20;
    Ms = [4, 8, 16, 32];
    color = ['r', 'g', 'b', 'y'];
    
    figure();
    hold on;
    scatter(SNRvalue,PerrorActuBinAp,'.','k');
    plot(SNRvalue,PerrorTheoBinAp,'k');
    for i = 1:length(Ms)
        % constellation for PSK
        theta = linspace(0, 2*pi, Ms(i)+1);
        theta = theta(1:Ms(i));
        constellationPSK = cos(theta) + 1j*sin(theta);
        
        % generate symbols for PSK
        symbolPSK = symbolGen(constellationPSK, nExample);
        
        % initializing theoretical and actual PSK Perror vector
        PerrorTheoPSK = zeros(1, length(SNRvalue));
        PerrorActuPSK = zeros(1, length(SNRvalue));

        for j = 1:length(SNRvalue)
            % simulation
            [nfBinPSK,receivedBinPSK,sfPSK] = scale(constellationPSK,symbolPSK,N0,j,false);

            % actual error
            PerrorActuPSK(j) = calcError(nfBinPSK, receivedBinPSK);

            % theoretical error - page 416 from the book
            PerrorTheoPSK(j) = 2*qfunc(sfPSK * sqrt(4/N0) * sin(pi/Ms(i)));
            
            
        end
        % plotting
        scatter(SNRvalue,PerrorActuPSK,'.',color(i));
        plot(SNRvalue,PerrorTheoPSK,color(i));
    end
    hold off;
    title('PSK bit error rate vs. SNR');
    ylabel('P_{error}');
    xlabel('SNR/bit (dB)');
    legend(["M=2","M=2 theoretical","M=4", "M=4 theoretical", "M=8", "M=8 theoretical", "M=16", "M=16 theoretical", "M=32", "M=32 theoretical"]);
    fprintf('Time taken for PSK with N = %d: %.2f s\n', nExample, toc());
    
%% Question 7
%DPSK
    tic;
    SNRvalue = -4:20;
    Ms = [4, 8, 16, 32];
    color = ['r', 'g', 'b', 'y'];
    
    figure();
    hold on;
    for i = 1:length(Ms)
        % constellation for DPSK
        theta = linspace(0, 2*pi, Ms(i)+1);
        theta = theta(1:Ms(i));
        constellationDPSK = cos(theta) + 1j*sin(theta);
        
        % generate symbols for DPSK
        symbolDPSK = symbolGen(constellationDPSK, nExample);
        
        % initializing theoretical and actual PSK Perror vector
        PerrorTheoDPSK = zeros(1, length(SNRvalue));
        PerrorActuDPSK = zeros(1, length(SNRvalue));

        for j = 1:length(SNRvalue)
            % simulation
            [nfBinDPSK,receivedBinDPSK,sfDPSK] = scale(constellationDPSK,symbolDPSK,N0,j,true);

            % actual error
            PerrorActuDPSK(j) = calcError(nfBinDPSK, receivedBinDPSK);

            % theoretical error - page 416 from the book
            PerrorTheoDPSK(j) = 2*qfunc(sfDPSK * sqrt(4/N0) * sin(pi/Ms(i)));
            
            
        end
        % plotting
        scatter(SNRvalue,PerrorActuDPSK,'.',color(i));
        plot(SNRvalue,PerrorTheoDPSK,color(i));
    end
    hold off;
    title('DPSK Perror vs. SNR');
    ylabel('P_{error}');
    xlabel('SNR/bit (dB)');
    legend(["M=4", "M=4 theoretical", "M=8", "M=8 theoretical", "M=16", "M=16 theoretical", "M=32 empirical", "M=32 theoretical"]);
    fprintf('Time taken for DPSK with N = %d: %.2f s\n', nExample, toc());
    
%% Question 8
%QAM
    tic;
    SNRvalue = -4:20;
    Ms = [4, 16, 32, 64];
    color = ['r', 'g', 'b', 'y'];
    
    figure();
    hold on;
    for i = 1:length(Ms)
        % constellation for QAM
        
        % square
        %if mod(log2(Ms(i)), 2) == 0
        if mod(log2(2), 2) == 0
            x = log2(Ms(i))-1
            y = log2(Ms(i))-1;
        
        % rectangle
        else
            x = log2(Ms(i))-2;
            y = Ms(i)/(x+1)-1;
        end
        xVec = -x/2:1:x/2;
        yVec = -y/2:1:y/2;

        [X,Y] = meshgrid(xVec,yVec);
        Y = Y * 1j;

        constellationQAM = X + Y;
        constellationQAM = reshape(constellationQAM,1,[]);
        
        % generate symbols for DPSK
        symbolQAM = symbolGen(constellationQAM, nExample);
        
        % initializing theoretical and actual PSK Perror vector
        PerrorTheoQAM = zeros(1, length(SNRvalue));
        PerrorActuQAM = zeros(1, length(SNRvalue));

        for j = 1:length(SNRvalue)
            % simulation
            [nfBinQAM,receivedBinQAM,sfQAM] = scale(constellationQAM,symbolQAM,N0,j,false);

            % energy
            Eng = mean(abs(nfBinQAM).^2);
            
            % actual error
            PerrorActuQAM(j) = calcError(nfBinQAM, receivedBinQAM);

            % theoretical error
            if mod(log2(Ms(i)), 2) == 0
                PerrorTheoQAM(j) = 1-(1-2*(1-1/sqrt(Ms(i)))*qfunc(sqrt(2*3*Eng/((Ms(i)-1) * N0))))^2;
            else
                PerrorTheoQAM(j) = 1-(1-2*qfunc(sqrt(2*3*Eng/((Ms(i)-1) * N0))))^2;
            end
            
        end
        % plotting
        scatter(SNRvalue,PerrorActuQAM,'.',color(i));
        plot(SNRvalue,PerrorTheoQAM,color(i));
    end
    hold off;
    title('QAM Perror vs. SNR');
    ylabel('P_{error}');
    xlabel('SNR/bit (dB)');
    legend(["M=4", "M=4 theoretical", "M=16", "M=16 theoretical", "M=32", "M=32 theoretical", "M=64", "M=64 theoretical"]);
    fprintf('Time taken for QAM with N = %d: %.2f s\n', nExample, toc());
%% Functions
function output = symbolGen(constellation, num)
    output = constellation(randsample(length(constellation),num,'true'));
end

function output = LSestimation(constellation, symbolNoisy)
    % distance from each symbol to every constallation
    distance = abs(constellation - symbolNoisy.');
    
    % index number of min value in each row
    [minDistance, iMin] = min(distance,[],2);
    
    % return the closest constellation
    output = constellation(iMin);
end

function [nf, received, scaleFactor] = scale(constellation, symbol, N0, SNRpBit, differential)
    
    M = length(constellation);
    N = size(symbol, 2);

    % desired average symbol energy
    % (from the SNR per bit, the noise power 
    %  and the number of symbols in the constellation)
    desiredAvgEnergy = 10^(SNRpBit/20) * N0/2 * log2(M);
    
    % average energy of constellation
    avgEnergy = mean(abs(constellation).^2);
    
    % scaling
    scaleFactor = sqrt(desiredAvgEnergy / avgEnergy);
    
    % scale base constellation to true constellation
    constellationScaled = constellation .* scaleFactor;
    
    % scale symbol to true transmitted scaled symbol
    symbolScaled = symbol .* scaleFactor;
    
    % generate AWGN
    %    (variance = N0/2)
    noise = sqrt(N0/2)/2 * (randn(1,N) + 1j*randn(1,N)); 
    
    % adding noise to symbol to get noisy symbol
    noisy = symbolScaled + noise;
    
    % differential scheme
    if differential == true
        %fprintf("in true\n");
        
        % phase of scaled symbol
        symbolScaledAng = atan2(imag(symbolScaled), real(symbolScaled));
        
        % return phase differential of noise-free scaled symbol
        nf = diff(symbolScaledAng);
        
        % phase of constellation
        constellationAng = atan2(imag(constellationScaled),real(constellationScaled));
        constellationAng = [constellationAng constellationAng+2*pi];
        constellationAngle = ([constellationAng-min(constellationAng) 4*pi])-2*pi; %min phase 0

        % phase of noisy symbol
        noisyAng = atan2(imag(noisy),real(noisy));
        
        % differential noisy
        noisyAng = diff(noisyAng);
        
        % call LSestimation and return estimated received vector
        received = LSestimation(constellationAngle, noisyAng);
    
    % non-differential scheme
    elseif differential == false
        %fprintf("in false\n");
        
        % return noise-free scaled symbol
        nf = symbolScaled;
        
        % call LSestimation and return estimated received vector
        received = LSestimation(constellationScaled, noisy);
    end

end

function output = calcError(estSym, trueSym)
    N = size(estSym, 2);
    threshold = 0.001;
    output = nnz(abs(estSym-trueSym) > threshold)/N;
end
