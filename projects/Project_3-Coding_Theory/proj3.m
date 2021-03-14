%% Error Control Coding Simulation Project
%   Bonny Wang Allister Liu Amy Leong
    clc
    clear
    close all
%% Q1
    %Use random matrix to test Galois field 2 Multiplication
    test1 = [1 1; 1 0];
    test2 = [1 0; 1 0];

    result_Multp = f2matrixMultp(test1,test2);
%% Q2
    %Use random matrix to test Galois field 2 Addition
    result_Add = f2matrixAdd(test1,test2);
%% Q3
    %Test corrupting str function
    cstr = corruptString('Hello World',0.5);
%% Rest Part Initalization
    G1 = [1,0,0,0,0,1,1,1;
          0,1,0,0,1,1,1,0;
          0,0,1,0,1,0,1,1;
          0,0,0,1,1,1,1,1];
    G2 = [1,0,0,0,1,1,1,1,0,1,1,0;
          0,1,0,0,1,0,0,1,1,1,1,0;
          0,0,1,0,1,1,0,1,1,0,1,1;
          0,0,0,1,1,0,1,0,1,1,1,1];
    k = 4;
%% Q4
    info = de2bi(0:16-1);
    info = flip(info,2);
  
    codeWord1 = f2matrixMultp(info,G1);
    codeWord2 = f2matrixMultp(info,G2);
    
    weight1 = sum(codeWord1,2);
    weight2 = sum(codeWord2,2);
    
    %dmin is equal to min weight
    dmin1 = min(weight1(weight1 > 0));
    dmin2 = min(weight2(weight2 > 0));
    
    maxError1 = (dmin1-1)/2;
    maxError2 = (dmin2-1)/2;
 %% Q5
    P1 = G1(:,k+1:end);
    P2 = G2(:,k+1:end);
    
    H1 = [ P1' eye(8-k)];
    H2 = [ P2' eye(12-k)];
    
    parity_Check1 = sum(f2matrixMultp(G1, H1'),'all')
    parity_Check2 = sum(f2matrixMultp(G2, H2'),'all')
    
    %They both equal to 0
    
%% Q6
    %Since there is only one error, same as an idnetity matrix
    all_Correctable1 = [zeros(1,8);eye(8)];
    
    all_Correctable2 = [zeros(1,12);eye(12)];
    for i = flip(1:11)
        %Add zeros before and identiy matrix after each position of 1 
        all_Correctable2 = [all_Correctable2; zeros(i,11-i) ones(i,1) eye(i)];
    end
        
    syndrome1 = f2matrixMultp(all_Correctable1,H1');
    syndrome2 = f2matrixMultp(all_Correctable2,H2');
    
    difftoReal1 = 2^(8-k) - size(syndrome1,1);
    % 7 smaller
    difftoReal2 = 2^(12-k) - size(syndrome2,1);
    % 177 smaller
    
%% Q7

    %Test with G2 with 2 message
    m = codeWord2(randsample(size(codeWord2, 1), 2, 'true'), : ); 
     
    erros = all_Correctable2(randsample(size(all_Correctable2, 1),2, 'true'), : );
     
    receivedM = f2matrixAdd(m,erros);    
    
    corrected = correctWords(all_Correctable2,syndrome2,receivedM, H2);
    % Original message is retrieved

%% Q8
    n = 100000;
    probs = linspace(0.001, 0.1, 20);
    
    m = info(randsample( size(info, 1), n, 'true' ), :); 
    coded_m1 = f2matrixMultp(m, G1 );
    coded_m2 = f2matrixMultp(m, G2 );
    
    uncodedRate = zeros(1,20);
    codedmRate1 = zeros(1,20);
    codedmRate2 = zeros( 1,20);
    
    for i = 1:20
        
        corruptedstr = corruptString(m, probs(i) );
        corrupted1 = corruptString(coded_m1, probs(i) );
        corrupted2 = corruptString(coded_m2, probs(i) );

        corrected1 = correctWords(all_Correctable1, syndrome1, corrupted1, H1 );
        corrected2 = correctWords(all_Correctable2, syndrome2, corrupted2, H2 );

        uncodedRate(i) = calcErroRate( m, corruptedstr );
        codedRate1(i) = calcErroRate( coded_m1, corrected1);
        codedRate2(i) = calcErroRate( coded_m2, corrected2);

    end
    
    figure();
    hold on;
    plot(probs, uncodedRate, 'DisplayName', 'UncodedM' );
    plot(probs, codedRate1, 'DisplayName', 'Coded G1' );
    plot(probs, codedRate2, 'DisplayName', 'Coded G2' );
    legend();
    title( "Plot of Bit Error Rate" );
    xlabel( "Probability" );
    ylabel( "Bit Error Rate" );
    hold off;

%% Functions

function result = f2matrixMultp(m1,m2)
    temp = m1*m2;
    
    %Since multiplication already followed the rule, just need to deal with
    %addition during matrix multiplication
    result = mod(temp,2);
end

function result = f2matrixAdd(m1,m2)
    temp = m1+m2;
    
    %Similar to binary
    result = mod(temp,2);
end

function output = corruptString(str,pError)
    [l,w] = size(str);    
    error = (rand( l,w ) <= pError);
    output = f2matrixAdd(str, error);
end

function corrected_Words = correctWords(errorWords, tsyndrom, receivedWords,H)
    
    % locate the error by multiply by H and compare with syndrom
    [ ~, index ] = ismember( f2matrixMultp(receivedWords,H'), tsyndrom, 'rows' ); 
    
    % Fill the last row will NaNs
    % If Index is 0, no correctable error, NaN will be the result
    index( index == 0 ) = size( errorWords,1)+1;
    errorWords = [errorWords; NaN( 1, size( errorWords, 2))];
    
    % 1 + 1 = 0 to eliminate the error
    corrected_Words = f2matrixAdd(errorWords(index,: ), receivedWords);
end

function rate = calcErroRate( m, predicted )
    m = m(:, 1:4);
    predicted = predicted(:, 1:4);

    bitErros = abs(m - predicted );
    
    % Uncorrectable count as a full word of errors.
    bitErros( isnan(bitErros)) = 1;
    
    rate = sum(bitErros, 'all' ) /numel(m);
end