% ECE408: Samuel Maltz
% CDMA Decoding Assignment
% Decodes a code-division multiple access (CDMA) message by filtering it,
% unspreading it, extracting it from the proper Walsh channel, demodulating
% it and decoding it from ASCII. To unspread the message, the correct
% pseudo-random (PN) sequence has to be used. To determine the correct
% shift of the PN sequence, a test is done on all possible shifts using the
% pilot-only frames.
clear; close all; clc;

load("Rcvd_Maltz.mat");    % received message Rcvd

H = hadamard(8);    % matrix of Walsh/Hadamard codes

pn = comm.PNSequence("Polynomial",[8 7 6 1 0],"InitialConditions",1, ...
    "Mask",[1 0 0 0 0 0 0 0],"SamplesPerFrame",255);    % PN sequence

% Root raised cosine filter.
h = [0.0038 0.0052 -0.0044 -0.0121 -0.0023 0.0143 0.0044 -0.0385 ...
    -0.0563 0.0363 0.2554 0.4968 0.6025 0.4968 0.2554 0.0363 -0.0563 ...
    -0.0385 0.0044 0.0143 -0.0023 -0.0121 -0.0044 0.0052 0.0038];

Nascii = 8;
Nsamp = 4;
Nwalsh = 8;
Nseq = pn.SamplesPerFrame;
Nfilt = length(h);
Nframe = Nsamp*Nseq;

% Walsh channels.
pilotchan = 0;
datachan = 5;

% Simulated pilot-only frame at the receiver.
pilot = zeros(ceil(Nseq/Nwalsh),1);
pilotbpsk = pskmod(pilot, 2);
pilotwalsh = pilotbpsk * H(pilotchan+1,:);

% Received pilot-only first frame.
Rcvdpilotunfilt = conv(h,Rcvd(1:Nframe-1));
Rcvdpilotdownsamp = Rcvdpilotunfilt(Nfilt:Nsamp:end);
p = zeros(Nseq,1);

% Loops through all possible shifts and calculates test metric.
for i = 1:Nseq
    Rcvdpilotspread = Rcvdpilotdownsamp .* pskmod(circshift(pn(),i),2).';
    p(i) = dot(pilotwalsh(1:end-1),Rcvdpilotspread);
end
    
figure;
plot(abs(p));
xlabel("n_{shift}");
ylabel("P(n_{shift})");
xlim([1 Nseq]);

% Index corresponding to max value of test metric is taken as proper PN
% sequence shift.
[~,shift] = max(abs(p));

% Decoding of CDMA message.
Rcvdunfilt = conv(h,Rcvd);                       % filtering
Rcvddownsamp = Rcvdunfilt(Nfilt:Nsamp:end);      % downsampling
Rcvdframes = reshape(Rcvddownsamp,Nseq,1,[]);

% Unspreading using the PN sequence.
Rcvdunspread = Rcvdframes .* pskmod(circshift(pn(),shift),2);

% Extracting the chips with encoded character.
Rcvddata = reshape(Rcvdunspread(1:Nascii*Nwalsh* ...
    floor(Nseq/(Nascii*Nwalsh)),1,:),Nwalsh,[],size(Rcvdunspread,3));

% Extracting the pilot from Walsh channel 0.
Rcvdpilot = pagemtimes(H(pilotchan+1,:),Rcvddata) / Nwalsh;

% Extracting the data from Walsh channel 5.
Rcvdunwalsh = squeeze(pagemtimes(H(datachan+1,:), ...
    (Rcvddata ./ mean(Rcvdpilot,2))) / Nwalsh);
Rcvdunwalsh = reshape(Rcvdunwalsh,Nascii,[]);

% Zeroing out low magnitude symbols.
Rcvdunwalsh(abs(Rcvdunwalsh) < 0.1) = 0;
Rcvddemod = pskdemod(Rcvdunwalsh,2);            % demodulation
msg = char(bit2int(Rcvddemod,Nascii,false));    % ASCII decoding
msg(msg == 0) = []  %#ok