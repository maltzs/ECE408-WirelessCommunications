% ECE408: Samuel Maltz
% OFDM-MIMO Assignment: Part 1
% Simulates flat fading 2x2 multiple-input multiple-output (MIMO) channels.
% Uses three different channels and three different equalization
% techniques: precoding, zero forcing and minimum mean squared error
% (MMSE). Attempts to maximize bit rate for a bit error rate threshold of
% 0.05.
clear; close all; clc;

Niter = 100;

fs = 1/4e-6;    % assumes a symbol period of 4us
Ntx = 2;
Nrv = 2;
snr = 20;

% Rows correspond to channels and columns correspond to equalization
% techniques.
M = [4 4 4; 4 4 2; 4 2 2];    % modulation order
N = fs*log2(M)*Ntx;           % number of bits in 1s

lambda = 1e-3;           % for MMSE equalization
fm = [1; 10; 100];       % max Doppler shifts of channels
ber = zeros(size(M));

for i = 1:length(fm)
    for j = 1:Niter
        % Creates 4 Rayleigh flat fading channels.
        chan = rayleigh(fm(i),4,Nrv,Ntx);
        H = chan(:,:,1);

        [U, ~, V] = svd(H);    % SVD for precoding equalization
        for k = 1:size(M,2)
            datatx = randi([0 1],N(i,k),1);
            symtx = pskmod(datatx,M(i,k),"InputType","bit");
			
			% Allocating alternating symbols for the two transmitters.
            x = reshape(symtx,Ntx,[]);
            
            % Multiplication on left by right singular vectors for
            % precoding equalization.
            if k == 1
                x = V * x;
            end

            % Transmission through the Rayleigh and additive white Gaussian
            % noise channel.
            y = H * x;
            for m = 1:Nrv
                y(m,:) = awgn(y(m,:),snr,'measured');
            end

            % Selection of equalization techniques.
            switch k
                case 1
                    symrv = U' * y;    % precoding
                case 2
                    symrv = pinv(H) * y;    % zero forcing
                case 3
                    symrv = H' * ...
                        (H * H' + lambda*eye(Nrv))^-1 * y;    % MMSE
            end
        
            datarv = pskdemod(symrv(:),M(i,k),"OutputType","bit");

            ber(i,k,j) = sum(datatx ~= datarv)/N(i,k);
        end
    end
end

table(fm,N(:,1),N(:,2),N(:,3),'VariableNames', ...
    ["Max Doppler shift (Hz)", "Precoding", "Zero forcing", "MMSE"])

table(fm,mean(ber(:,1,:),3),mean(ber(:,2,:),3),mean(ber(:,3,:),3), ...
    'VariableNames',["Max Doppler shift (Hz)", "Precoding", ...
    "Zero forcing", "MMSE"])