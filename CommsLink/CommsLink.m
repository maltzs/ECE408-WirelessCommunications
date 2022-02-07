% ECE408: Samuel Maltz
% Communications Link Assignment
% A BER script for a wireless link simulation. This script performs 3 tasks
% which are taken from the 2 parts of Prof. Keene's final project of
% ECE300.
% Case 1: BPSK, 4QAM and 16QAM on no ISI channel for all SNR values (Part 1)
% Case 2: BPSK on moderate and severe ISI channels with equalization (Part 1)
% Case 3: BPSK on moderate and severe ISI channels with equalization and coding (Part 2)
clear; close all; clc;
% For the final version of this project, you must use these 3
% parameter. You will likely want to set numIter to 1 while you debug your
% link, and then increase it to get an average BER.
numIter = 1000;  % The number of iterations of the simulation
nSym = 1000;    % The number of symbols per packet
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);

Marys = {[2; 4; 16], 2, 2};        % The M-ary number, 2 corresponds to binary modulation, different M-ary numbers used based on task

% 3 types of channels including: no channel, somewhat invertible channel
% impulse response (moderate ISI) and not so invertible (severe ISI).
chan = {1, [1 0.2 0.4], [0.227 0.460 0.688 0.460 0.227]};
chanTasks = {1, [2 3], [2 3]};  % Different channels used based on task
chanStrings = ["No ISI", "Moderate ISI", "Severe ISI"];

nTrain = 100;  % Number of training symbols for equalization

% Parameters for BCH code used in task 3 (4/7 code)
m = 3;
n = 2^m-1;
k = 4;

% Time-varying Rayleigh multipath channel, try it if you dare. Or take
% wireless comms next semester.
% ts = 1/1000;
% chan = rayleighchan(ts,1);
% chan.pathDelays = [0 ts 2*ts];
% chan.AvgPathGaindB = [0 5 10];
% chan.StoreHistory = 1; % Uncomment if you want to be able to do plot(chan)

for task = 1:3
    M = Marys{task};

    % Create a vector to store the BER computed during each iteration
    berVec = zeros(numIter, lenSNR, length(chanTasks{task}));

    % Create a vector to store the bit rate computed for each M-ary number
    br = zeros(length(M),1);
    
    for q = 1:length(chanTasks{task})
        for p = 1:length(M)
            % Run the simulation numIter amount of times
            for i = 1:numIter
                nBits = nSym*log2(M(p));

                % Reduce numBits by extra code bits for task 3
                if task == 3
                    nBits = floor(nBits/n)*k;
                end

                bits = randi([0 1], 1, nBits);     % Generate random bits
                % New bits must be generated at every
                % iteration
            
                % If you increase the M-ary number, as you most likely will, you'll need to
                % convert the bits to integers. See the BIN2DE function
                % For binary, our MSG signal is simply the bits
                % For task 3, encode bits using BCH encoding
                if task == 3
                    msg = gf(reshape(bits,[],k));
                    msg = bchenc(msg,n,k);
                    msg = double(msg.x(:));
                else
                    msg = bits';
                end
            
                for j = 1:lenSNR % one iteration of the simulation at each SNR Value
                    % BPSK if M == 2, QAM otherwise
                    if M(p) == 2
                        tx = pskmod(msg,M(p));
                    else
                        tx = qammod(msg,M(p),'InputType','bit');  % BPSK modulate the signal
                    end
            
                    if isequal(chan{chanTasks{task}(q)},1)
                        txChan = tx;
                    elseif isa(chan{chanTasks{task}(q)},'channel.rayleigh')
                        reset(chan) % Draw a different channel each iteration
                        txChan = filter(chan{chanTasks{task}(q)},tx);
                    else
                        txChan = filter(chan{chanTasks{task}(q)},1,tx);  % Apply the channel
                    end
             
                    txNoisy = awgn(txChan,SNR_Vec(j),'measured'); % Add AWGN
            
                    % Equalization for tasks 2 and 3
                    if task ~= 1
                        equalizer = comm.LinearEqualizer("Algorithm","RLS",...
                        "NumTaps",3,"Constellation",pskmod(0:M(p)-1,M(p)), ...
                        "ReferenceTap",1,"InitialInverseCorrelationMatrix",1);
        
                        [txEq,err] = equalizer(txNoisy,tx(1:nTrain));

                        % Scatter plot of equalization
                        if task == 2 && q == 1 && i == 1 && SNR_Vec(j) == 12
                            c = [1 1 0; 0 0 1];
                            cons = pskmod([0 M-1],M);
                            figure;
                            scatter([real(txNoisy(nTrain:end)) real(txEq(nTrain:end))], ...
                                [imag(txNoisy(nTrain:end)) imag(txEq(nTrain:end))],20,c,'.');
                            hold on;
                            scatter(real(cons), imag(cons),40,'r','+');
                            hold off;
                            axis equal;
                            grid on;
                            set(gca,'Color','k');
                            set(gca,'GridColor','w');
                            xlabel("In-phase Amplitude");
                            ylabel("Quadrature Amplitude");
                            legend(["Before Equalization", "After Equalization"],'TextColor','w');
                        end
                    else
                        txEq = txNoisy;
                    end
            
                    if M(p) == 2
                        rx = pskdemod(txEq,M(p));
                    else
                        rx = qamdemod(txEq,M(p),'OutputType','bit'); % Demodulate
                    end
            
                    % For task 3, decode bits using BCH decoder
                    if task == 3
                        dec = gf(reshape(rx,[],n));
                        dec = bchdec(dec,n,k);
                        rxMSG = double(dec.x(:))';
                    else
                        rxMSG = rx';
                    end
                    
                    % Compute and store the BER for this iteration
                    [~,berVec(i,j,q)] = biterr(bits, rxMSG);  % We're interested in the BER, which is the 2nd output of BITERR
            
                end  % End SNR iteration
            end      % End numIter iteration        
        
            % Compute and plot the mean BER
            ber = mean(berVec,1);
    
            % Compute number of useable bits
            if task == 1
                % Task 1: All bits
                nUseablebits = nBits;
            elseif task == 2
                % Task 2: All bits minus training bits
                nUseablebits = nBits-log2(M(p))*nTrain;
            else
                % Task 3: All bits minus the training bits that are message
                % bits (nBits already does not contain extra code bits)
                nUseablebits = nBits-floor(log2(M(p))*nTrain/n)*k;
            end
    
            br(p) = nUseablebits/nSym;
    
            % Results: BER curves for task 1, BER and bit rate at 12 dB SNR on
            % moderate ISI and severe ISI channels for tasks 2 and 3
            figure;
            semilogy(SNR_Vec, ber(:,:,q));
            
            % Compute the theoretical BER for this scenario
            % THIS IS ONLY VALID FOR BPSK!
            % YOU NEED TO CHANGE THE CALL TO BERAWGN FOR DIFF MOD TYPES
            % Also note - there is no theoretical BER when you have a multipath channel
            % Use BPSK if M == 2, QAM otherwise
            if M(p) == 2
                berTheory = berawgn(SNR_Vec,'psk',2,'nondiff');
            else
                berTheory = berawgn(SNR_Vec-10*log10(log2(M(p))),'qam',M(p)); % berawgn accepts EbNo, not SNR
            end
            hold on
            semilogy(SNR_Vec,berTheory,'r');
            xlabel("SNR (dB)");
            ylabel("BER");
            legend('BER', 'Theoretical BER', 'location', 'southwest');
        end          % End M-ary iteration
    
        % Summary of bit rates vs. M-ary number
        table(repmat(chanStrings(chanTasks{task}(q)),length(M),1),M,br,'VariableNames',["Channel", "M", "Bit Rate"])
    end              % End chan iteration
end                  % End task iteration