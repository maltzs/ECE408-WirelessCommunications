% ECE408: Samuel Maltz
% Alamouti Codes Assignment
% Recreates Figure 4 in "A Simple Transmit Diversity Technique
% for Wireless Communications" by Alamouti by simulating a Rayleigh channel
% and the three schemes used.
% Case 1: No diversity (1 transmitter, 1 receiver)
% Case 2: Maximal-ratio receiver combining (MRRC) with subcases (1
% transmitter, 2 receivers) and (1 transmitter, 4 receivers)
% Case 3: New scheme with subcases (2 transmitters, 1 receiver) and (2
% transmitters, 2 receivers)
clear; close all; clc;

nIter = 1e3;    % Iterations
nSym = 2^10;    % Message length
SNR = 0:2:48;

% Number of transmitters and receivers for the three cases.
nTx = {1, [1 1], [2 2]};
nRx = {1, [2 4], [1 2]};

M = 2;  % BPSK

% 900MHz carrier frequency, Doppler shift caused by 60mi/hr driving car.
fc = 900e6;
v = 26.66;
c = 3e8;
fm = fc*v/c;

BER = cell(1,3);
for i = 1:length(nTx)
    BER{i} = zeros(nIter,length(SNR),length(nTx{i}));
end

% Iterations loop
for i = 1:nIter
    bits = randi([0 1],nSym,1);
    tx = pskmod(bits,M);

    % Cases loop
    for j = 1:length(nTx)
        % Subcases loop
        for k = 1:length(nTx{j})
            % If running new technique, converts original symbols into 2
            % columns, one for each transmitter and computes appropriate 
            % sent symbols.
            if j == 3
                tx1 = reshape(tx,nTx{j}(k),[]);
                tx1 = [reshape([tx1(1,:); -conj(tx1(2,:))],nSym,1) ...
                    reshape([tx1(2,:); conj(tx1(1,:))],[],1)];
            else
                tx1 = tx;
            end

            % Creates Rayleigh channels, assumes same channel response for
            % consecuitive symbols.
            rayChan1 = rayleigh(fc, fm, nSym/nTx{j}(k), ...
                nTx{j}(k)*nRx{j}(k));
            rayChan = repelem(rayChan1,nTx{j}(k),1,1);
            
            % Applies channels.
            txRayChan = rayChan.*repmat(tx1,1,nRx{j}(k));

            % SNR loop
            for m = 1:length(SNR)
                % Applies AWGN.
                r = awgn(txRayChan,SNR(m),"measured");  
                
                % If running new scheme, combine received symbols as per
                % scheme.
                if j == 3
                    % Tried to only use one branch regardless of number of
                    % receivers, could not come up with a general formula.
                    if k == 1
                        % Sum received symbols as per scheme, reshape r
                        % into nSym/2 x 2 matrix with rows i being r0 r1
                        % for the ith symbol pair.
                        r = reshape(sum(r,2),nTx{j}(k),[]).';

                        % Compute stilde0s per scheme.
                        stilde0 = conj(rayChan1).*r;
                        stilde0 = sum([stilde0(:,1) conj(stilde0(:,2))],2);
    
                        % Compute stilde1s per scheme.
                        stilde1 = conj(rayChan1).*[r(:,2) r(:,1)];
                        stilde1 = sum([-conj(stilde1(:,1)) ...
                            stilde1(:,2)],2);
                    else
                        % Sum received symbols as per scheme, reshape r
                        % into nSym/2 x 4 matrix with rows i being r0 r1 r2
                        % r3 for the ith symbol pair.
                        r = reshape(sum(reshape(r.',nTx{j}(k),[])), ...
                            nTx{j}(k)*nRx{j}(k),[]).';
                        r(:,2:3) = fliplr(r(:,2:3));

                        % Computes stilde0s per scheme.
                        stilde0 = conj(rayChan1).*r;
                        stilde0 = sum([stilde0(:,1) conj(stilde0(:,2)) ...
                            stilde0(:,3) conj(stilde0(:,4))],2);
    
                        % Computes stilde1s per scheme.
                        stilde1 = conj(rayChan1).*[r(:,2) r(:,1) ...
                            r(:,4) r(:,3)];
                        stilde1 = sum([-conj(stilde1(:,1)) stilde1(:,2) ...
                            -conj(stilde1(:,3)) stilde1(:,4)],2);
                    end

                    % Interleaves the stilde0s and stilde1s.
                    stilde = reshape([stilde0 stilde1].',[],1);
                else
                    % stilde without new scheme.
                    stilde = sum(conj(rayChan).*r,2);
                end

                rx = pskdemod(stilde,M);
    
                [~, BER{j}(i,m,k)] = biterr(bits,rx);
            end
        end
    end
end

markers = {{'o'}, {'v';'s'}, {'d';'^'}};
legends = strings(1,sum(cellfun(@numel, markers)));
k = 1;

figure;
for i = 1:length(nTx)
    BERavg = reshape(mean(BER{i},1),length(SNR),length(nTx{i}));
    h = semilogy(SNR,BERavg);
    set(h,{'Marker'},markers{i});
    for j = 1:length(nTx{i})
        switch i
            case 1
                legends(k) = "no diversity";
            case 2
                legends(k) = "MRRC";
            case 3
                legends(k) = "new scheme";
        end

        legends(k) = legends(k) + " (" + nTx{i}(j) + " Tx, " + ...
            nRx{i}(j) + " Rx)";
        k = k+1;
    end
    hold on;
end
hold off;
grid on;
xlim([0 50]);
ylim([1e-6 1e0]);
xticks(5:5:50);
xlabel("SNR (dB)");
ylabel("P_b, bit error rate (BER)");
legend(legends);


% Simulates n Rayleigh channels of length N assuming carrier frequency of
% fc and max Doppler shift of fm. Algorithm from "Wireless Communications:
% Principles and Practice," 1st Edition by Rappaport.
function r = rayleigh(fc, fm, N, n)
    f = fc + linspace(-fm,fm,N)';

    g = randn(N/2,n) + 1j*randn(N/2,n);
    g = [conj(flipud(g)); g];
    
    % Fading spectrum, uses linearized values at S(1) and S(end) due to
    % infinity values.
    S = 1.5./(pi*fm*sqrt(1-((f(2:end-1)-fc)/fm).^2));
    dS = S(2)-S(1);
    S = [S(1)-dS; S; S(end)-dS];
    
    g = sqrt(S).*g;
    gi = ifft(real(g));
    gq = ifft(imag(g));
    
    % Rayleigh is magnitude of complex Gaussian random variable
    r = sqrt(abs(gi).^2+abs(gq).^2);
end