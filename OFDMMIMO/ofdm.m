% ECE408: Samuel Maltz
% OFDM-MIMO Assignment: Part 2
% Simulates transmission of 802.11a orthogonal frequency-division
% multiplexing (OFDM) data symbols through frequency-selective channels.
% Uses three single path channels and two equalization techniques: zero
% forcing and minimum mean squared error (MMSE). Attempts to maximize bit
% rate for a bit error rate threshold of 0.05.
clear; close all; clc;

Niter = 100;

Nsamp = 80;
Nsc = 64;
Ndatasc = 48;
Ts = 4e-6;
fs = Nsamp/Ts;
snr = 20;

% Rows correspond to channels and columns correspond to equalization
% techniques.
M = [64 64; 16 16; 2 2];    % modulation order
N = log2(M)*Ndatasc/Ts;     % number of bits in 1s

% Data and pilot subcarriers.
datasc = -26:26;
pilotsc = [-21 -7 7 21];
datasc(any(datasc == pilotsc') | datasc == 0) = [];

% Frequency-selective channels.
chan = {[0.04 -0.05 0.07 -0.21 -0.5 0.72 0.36 0 0.21 0.03 0.07], ...
    [0.407 0.815 0.407], [0.227 0.460 0.688 0.460 0.227]};
ber = zeros([size(M) Niter]);

for i = 1:length(chan)
    H = freqz(chan{i},1,Nsc,fs,'whole');    % frequency response of channel
    for j = 1:Niter
        for k = 1:size(M,2)
            datatx = randi([0 1],N(i,k),1);
            if M(i,k) <= 4
                symtx = pskmod(datatx,M(i,k),"InputType","bit");
            else
                symtx = qammod(datatx,M(i,k),"InputType","bit");
            end
            
            % Placement of data on data subcarriers.
            gridtx = zeros(Nsc,length(symtx)/Ndatasc);
            gridtx(datasc+Nsc/2+1,:) = reshape(symtx,length(datasc),[]);
            
            % Placement of pilot on pilot subcarriers.
            pn = comm.PNSequence("Polynomial",[7 4 0], ...
                "InitialConditions",1,"Mask",[0 0 0 0 0 0 1], ...
                "SamplesPerFrame",length(symtx)/Ndatasc);
            pilotseq = [1; 1; 1; -1] .* pskmod(circshift(pn(), ...
                length(symtx)/Ndatasc-7),2,"InputType","bit").';
            gridtx(pilotsc+Nsc/2+1,:) = pilotseq;
            
            % Inverse fast Fourier transform and cyclic prefix prepending.
            gridtx = ifft(ifftshift(gridtx,1));
            gridtx = [zeros(Nsamp-Nsc,length(symtx)/Ndatasc); ...
                gridtx];    %#ok
            gridtx(1:Nsamp-Nsc,:) = gridtx(Nsc+1:end,:);
            
            sigman2 = var(gridtx(:));    % for MMSE equalization

            % Transmission through frequency-selective and additive white
            % Gaussian noise channel.
            [gridrv, sigmaw2] = awgn(filter(chan{i},1,gridtx(:)),snr, ...
                'measured');    % sigmaw2 for MMSE equalization

            % Removal of cyclic prefix and fast Fourier transform.
            gridrv = reshape(gridrv,Nsamp,[]);
            gridrv(1:Nsamp-Nsc,:) = [];
            gridrv = fft(gridrv);

            % Selection of equalization technique.
            if k == 1
                y = gridrv .* (1 ./ H);    % zero forcing
            else
                y = gridrv .* (conj(H) ./ ...
                    (abs(H).^2+sigmaw2/sigman2));    % MMSE
            end

            % Extraction of data from data subcarriers. 
            y = fftshift(y,1);
            symrv = reshape(y(datasc+Nsc/2+1,:),[],1);
            if M(i,k) <= 4
                datarv = pskdemod(symrv,M(i,k),"OutputType","bit");
            else
                datarv = qamdemod(symrv,M(i,k),"OutputType","bit");
            end
    
            ber(i,k,j) = sum(datatx ~= datarv)/N(i,k);
        end 
    end
end

table((1:3)',N(:,1),N(:,2),'VariableNames', ...
    ["Channel", "Zero forcing", "MMSE"])

table((1:3)',mean(ber(:,1,:),3),mean(ber(:,2,:),3),'VariableNames', ...
    ["Channel", "Zero forcing", "MMSE"])