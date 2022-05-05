% ECE408: Samuel Maltz
% OFDM-MIMO Assignment: Part 3
% Simulates orthogonal frequency-division multiplexing multiple-input
% multiple-output (OFDM-MIMO) over a frequency-selective 2x2 MIMO channel.
% Demonstrates how OFDM causes the frequency-selective channel to be split
% into a collection of narrow flat channels on which MIMO equalization can
% be applied. Zero forcing equalization is used.
clear; close all; clc;

Niter = 100;

Ntx = 2;
Nrv = 2;

Nsamp = 80;
Nsc = 64;
Ndatasc = 48;
Ts = 4e-6;
fs = Nsamp/Ts;
snr = 20;

M = 4;                         % modulation order
N = log2(M)*Ntx*Ndatasc/Ts;    % number of bits in 1s

% Data and pilot subcarriers.
datasc = -26:26;
pilotsc = [-21 -7 7 21];
datasc(any(datasc == pilotsc') | datasc == 0) = [];

fm = 1;
ber = zeros(1,Niter);

for i = 1:Niter
    % Creates 4 frequency-selective Rayleigh channels.
    chan = rayleigh(fm,Nsamp-Nsc,Nrv,Ntx);
    
    % Computes pseudo inverse for all subcarrier MIMO channel frequency
    % response matrices to use in zero forcing equalization.
    Hpi = zeros(Nrv,Ntx,Nsc);
    for j = 1:Nrv
        for k = 1:Ntx
            Hpi(j,k,:) = freqz(squeeze(chan(j,k,:)),1,Nsc,fs,'whole');
        end
    end
    
    for j = 1:Nsc
        Hpi(:,:,j) = pinv(Hpi(:,:,j));
    end
    
    datatx = randi([0 1],N,1);
    if M <= 4
        symtx = pskmod(datatx,M,"InputType","bit");
    else
        symtx = qammod(datatx,M,"InputType","bit");
    end
    
    % Placement of data on data subcarriers.
    gridtx = zeros(Nsc,length(symtx)/Ndatasc);
    gridtx(datasc+Nsc/2+1,:) = reshape(symtx,length(datasc),[]);
    
    % Placement of pilot on pilot subcarriers.
    pn = comm.PNSequence("Polynomial",[7 4 0],"InitialConditions", ...
        1,"Mask",[0 0 0 0 0 0 1],"SamplesPerFrame",length(symtx)/Ndatasc);
    pilotseq = [1; 1; 1; -1] .* pskmod(circshift(pn(), ...
        length(symtx)/Ndatasc-7),2,"InputType","bit").';
    gridtx(pilotsc+Nsc/2+1,:) = pilotseq;
    
    % Inverse fast Fourier transform and cyclic prefix prepending.
    gridtx = ifft(ifftshift(gridtx,1));
    gridtx = [zeros(Nsamp-Nsc,length(symtx)/Ndatasc); gridtx];    %#ok
    gridtx(1:Nsamp-Nsc,:) = gridtx(Nsc+1:end,:);

    % Allocating alternating OFDM symbols for the two transmitters.
    gridtx = reshape(permute(reshape(gridtx,Nsamp,Ntx,[]),[1 3 2]),[],Ntx);
    
    % Transmission through the Rayleigh and additive white Gaussian noise
    % channel.
    gridrv = zeros(size(gridtx));
    for j = 1:Nrv
        for k = 1:Ntx
            gridrv(:,j) = gridrv(:,j) + filter(squeeze(chan(j,k,:)), ...
                1,gridtx(:,k));
        end
        gridrv(:,j) = awgn(gridrv(:,j),snr,'measured');
    end
    
    % Removal of cyclic prefix and fast Fourier transform.
    gridrv = reshape(gridrv,Nsamp,[]);
    gridrv(1:Nsamp-Nsc,:) = [];
    gridrv = fft(gridrv);

    % Rearranging of data and zero forcing equalization on each subcarrier
    % separately.
    gridrv = permute(reshape(reshape(gridrv,[],Nrv).',Nrv,Nsc,[]),[1 3 2]);
    y = pagemtimes(Hpi,gridrv);
    
    % Extraction of data from data subcarriers.
    y = fftshift(y,3);
    symrv = reshape(permute(y(:,:,datasc+Nsc/2+1),[3 1 2]),[],1);
    if M <= 4
        datarv = pskdemod(symrv,M,"OutputType","bit");
    else
        datarv = qamdemod(symrv,M,"OutputType","bit");
    end
    
    ber(i) = sum(datatx ~= datarv)/N;
end

table(N,mean(ber),'VariableNames',["Bit rate", "Bit error rate"])