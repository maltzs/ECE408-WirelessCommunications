% ECE408: Samuel Maltz
% LTE Downlink Simulation Project
% Simulates a Long-Term Evolution (LTE) downlink transmission for 8
% reference channels listed in 3rd Generation Partnership Project (3GPP)
% technical specification (TS) 36.101 Annex A.3.3.1 Tables A.3.3.1-1 to
% A.3.3.1-3:
% R.2: Quadrature phase shift keying (QPSK), 50 resource blocks (RBs)
% R.3: 16-Quadrature amplitude modulation (QAM), 50 RBs
% R.4: QPSK, 6 RBs
% R.5: 64QAM, 15 RBs
% R.6: 64QAM, 25 RBs
% R.7: 64QAM, 50 RBs
% R.8: 64QAM, 75 RBs
% Calculates bit rate (BR) and bit error rate (BER) curves of physical
% downlink shared channel (PDSCH) transmission. Therefore, this program
% only simulates channels and signals which are dependent on the PDSCH.
clear; close all; clc;

NFrames = 1000;    % number of frames
snr = -10:2:20;    % signal-to-noise ratio (SNR) values
rc = 2:8;          % reference channel numbers

BR = zeros(1,length(rc));
AvgBER = zeros(length(snr),length(rc));

% Reference channel loop.
for i = 1:length(rc)
    BER = zeros(length(snr),NFrames);

    enbref = lteRMCDL("R." + rc(i));    % reference channel
    
    enbref.NRC = rc(i);
    enbref.DCIFormat = enbref.PDSCH.DCIFormat;

    % enbref.PDSCH.NSoftbits is for user equipment category 2 as listed in
    % 3GPP TS 36.306 Section 4.1 Table 4.1-1.          
    enbref.PDSCH.NSoftbits = 1237248;
    enbref.PDSCH.DuplexMode = enbref.DuplexMode;
    
    % SNR loop.
    for j = 1:length(snr)
        % Frame loop.
        for NFrame = 0:NFrames-1
            frame = [];
            payloadframe = [];
            payloadframercv = [];

            % Error decoding the physical downlink control channel (PDCCH).
            pdccherr = 0;

            enb = enbref;    % resets the reference channel
        
            % Subframe loop.
            for NSubframe = 0:enb.TotSubframes-1
                enb.NFrame = NFrame;
                enb.NSubframe = NSubframe;
        
                payload = randi([0 1], ...
                    enb.PDSCH.TrBlkSizes(enb.NSubframe+1),1);
                payloadframe = [payloadframe; payload];          %#ok
        
                grid = lteResourceGrid(enb);
                
                [pdschind, pdschinfo] = ltePDSCHIndices(enb,enb.PDSCH, ...
                    enb.PDSCH.PRBSet);
                pdschsym = generatePDSCH(enb,pdschinfo.G,payload);
                
                % 36.211-6.3.5: Mapping PDSCH symbols to resource elements
                grid(pdschind) = pdschsym;
                
                pcfichind = ltePCFICHIndices(enb);
                pcfichsym = generatePCFICH(enb);
                
                % 36.211-6.7.4: Mapping physical control format indicator
                % channel (PCFICH) symbols to resource elements
                grid(pcfichind) = pcfichsym;
                
                pdcchind = ltePDCCHIndices(enb);
                pdcchinfo = ltePDCCHInfo(enb);
                pdcchsym = generatePDCCH(enb, ...
                    pdcchinfo.MTot);
                
                % 36.211-6.8.5: Mapping PDCCH symbols to resource elements
                grid(pdcchind) = pdcchsym;
                
                % 36.211-6.11.1.1: Primary synchronization signal (PSS)
                % generation
                pssind = ltePSSIndices(enb);
                psssym = ltePSS(enb);
                
                % 36.211-6.11.1.2: Mapping primary synchronization signal
                % to resource elements
                grid(pssind) = psssym;
                
                % 36.211-6.11.2.1: Secondary synchronization signal (SSS)
                % generation
                sssind = lteSSSIndices(enb);
                ssssym = lteSSS(enb);
                
                grid(sssind) = ssssym;
        
                frame = [frame grid];   %#ok
            end 
        
            % 36.211-6.12: Orthogonal frequency-division multiplexing
            % (OFDM) baseband signal generation
            [ofdm, info] = lteOFDMModulate(enb,frame);
    
            ofdm = awgn(ofdm,snr(j),"measured");
    
            % Calculation of timing offset using PSS and SSS. Should be 0.
            % If not, do not decode frame.
            timingoffset = lteDLFrameOffset(enb,ofdm);
            if timingoffset == 0   
                framercv = lteOFDMDemodulate(enb,ofdm(timingoffset+1:end));
                
                l = size(framercv,2)/enb.TotSubframes;
                for NSubframe = 0:enb.TotSubframes-1
                    enb.NSubframe = NSubframe;
            
                    grid = framercv(:,NSubframe*l+1:(NSubframe+1)*l);
                    
                    % Decode PCFICH to recover control format indicator
                    % (CFI).
                    pcfichind = ltePCFICHIndices(enb);
                    pcfichsym = grid(pcfichind);
                    enb.CFI = decodePCFICH(enb,pcfichsym);

                    % If CFI is not correct, stop decoding frame.
                    if enb.CFI ~= enbref.CFI
                        pdccherr = 1;
                        break
                    end
    
                    % Decode PDCCH to recover downlink control information
                    % (DCI).
                    pdcchind = ltePDCCHIndices(enb);
                    pdcchsym = grid(pdcchind);
                    dcircv = decodePDCCH(enb,pdcchsym);
            
                    % 36.213-7.1.6: Recover resource allocation of PDSCH
                    PRBSet = (0:enb.NDLRB-1)';
                    PRB = repelem(dcircv.Allocation.Bitmap-'0', ...
                        ceil(enb.NDLRB/length(dcircv.Allocation.Bitmap)));
                    PRB(enb.NDLRB+1:end) = [];
                    PRBSet(~PRB) = [];
                    pdschind = ltePDSCHIndices(enb,enb.PDSCH,PRBSet);
                    pdschsym = grid(pdschind);

                    % For subframe 5 which does not have PDSCH symbols.
                    if isempty(pdschsym)
                        payloadrcv = [];
                    else
                        payloadrcv = decodePDSCH(enb,dcircv,pdschsym);

                        % If there was an error recovering the DCI, stop
                        % decoding the frame.
                        if payloadrcv == -1
                            pdccherr = 1;
                            break
                        end
                    end
            
                    payloadframercv = ...
                        [payloadframercv; payloadrcv];    %#ok
                end
            end
    
            % If there was an error decoding the PDCCH which led to a too
            % short payload, append with zeros.
            if pdccherr == 1 || ...
                    length(payloadframercv) < length(payloadframe)
                payloadframercv(end+1:length(payloadframe),1) = ...
                    zeros(length(payloadframe)-length(payloadframercv),1);
            end
    
            % If there was an error decoding the PDCCH which led to a too
            % long payload, truncate.
            if length(payloadframercv) > length(payloadframe)
                payloadframercv(length(payloadframe)+1:end) = [];
            end
    
            BER(j,NFrame+1) = sum(payloadframercv ~= payloadframe) ...
                /length(payloadframe);
        end
    end
    
    BR(i) = length(payloadframe)/10e-3;    % 1 frame is 10ms    
    AvgBER(:,i) = mean(BER,2);
end

table("R." + rc',BR','VariableNames',["Channel","BR"])

figure;
semilogy(snr,AvgBER);
xlabel("SNR (dB)");
ylabel("Average BER");
legend("R." + rc,"Location","southwest");


% Generates PDSCH symbols according to 3GPP TS 36.211 and TS 36.212.
function precodeblk = generatePDSCH(enb,G,payloadblk)    
    % 36.212-5.3.2.1: Cyclic redundancy check (CRC) attachment
    crcblk = lteCRCEncode(payloadblk,"24A");
    
    % 36.212-5.3.2.2: Code block segmentation
    segblk = lteCodeBlockSegment(crcblk);
    
    % 36.212-5.3.2.3: Channel coding
    turboblk = lteTurboEncode(segblk);
    
    % 36.212-5.3.2.4 & 36.212-5.3.2.5: Rate matching and concatenation
    matchedblk = lteRateMatchTurbo(turboblk,G,enb.PDSCH.RV,enb.PDSCH);
    
    % 36.211-6.3.1: Scrambling
    pdschscrseq = ltePDSCHPRBS(enb,enb.PDSCH.RNTI,0,length(matchedblk));
    scrblk = xor(pdschscrseq,matchedblk);
    
    % 36.211-6.3.2: Modulation
    modblk = lteSymbolModulate(scrblk,enb.PDSCH.Modulation);
    
    % 36.211-6.3.3: Layer mapping
    layermapblk = lteLayerMap(modblk,enb.PDSCH.NLayers);
    
    % 36.211-6.3.4: Precoding
    precodeblk = lteDLPrecode(enb,enb.PDSCH,layermapblk);
end

% Generates PCFICH symbols according to 3GPP TS 36.211 and TS 36.212.
function precodeblk = generatePCFICH(enb)
    % 36.212-5.3.4.1: Channel coding
    cfiblk = lteCFI(enb);

    % 36.211-6.7.1: Scrambling
    pcfichscrseq = ltePCFICHPRBS(enb,length(cfiblk));
    scrblk = xor(pcfichscrseq,cfiblk);

    % 36.211-6.7.2: Modulation
    modblk = lteSymbolModulate(scrblk,"QPSK");

    % 36.211-6.7.3: Layer mapping and precoding
    layermapblk = lteLayerMap(modblk,enb.PDSCH.NLayers);
    precodeblk = lteDLPrecode(enb,enb.PDSCH,layermapblk);
end

% Generates PDCCH symbols according to 3GPP TS 36.211 and TS 36.212.
function multiplexblk = generatePDCCH(enb,MTot)
    % Computes modulation and coding scheme (MCS) index of 3GPP TS 36.213.
    switch enb.NRC
        case 2
            mcs = 5;
        case 3
            if enb.NSubframe == 0
                mcs = 14;
            else
                mcs = 15;
            end
        case 4
            if enb.NSubframe == 0
                mcs = 0;
            else
                mcs = 4;
            end
        case 5
            if enb.NSubframe == 0
                mcs = 21;
            else
                mcs = 25;
            end
        case 6
            if enb.NSubframe == 0
                mcs = 23;
            else
                mcs = 25;
            end
        case 7
            if enb.NSubframe == 0
                mcs = 25;
            else
                mcs = 26;
            end
        case 8
            if enb.NSubframe == 0
                mcs = 26;
            else
                mcs = 27;
            end
    end
    
    % For resource allocation of PDSCH from 3GPP TS 36.213.
    if enb.NDLRB <= 10
        P = 1;
    elseif enb.NDLRB <= 26
        P = 2;
    elseif enb.NDLRB <= 63
        P = 3;
    else
        P = 4;
    end
    
    % 36.212-5.3.3.1: DCI generation
    dciin.DCIFormat = enb.DCIFormat;
    dciin.AllocationType = 0;

    % All resource blocks are used for PDSCH except in subframe 5 when none
    % are.
    dciin.Allocation.Bitmap = dec2bin(2^ceil(enb.NDLRB/P)-1);
    if enb.NSubframe == 5
        dciin.Allocation.Bitmap = dec2bin(0,ceil(enb.NDLRB/P));
    end
    dciin.ModCoding = mcs;
    dciin.RV = enb.PDSCH.RV;
    dciin.TPCPUCCH = 0;
    [~, dciblk] = lteDCI(enb,dciin);
    
    % 36.212-5.3.3.2: CRC attachment
    crcblk = lteCRCEncode(dciblk,"16",enb.PDSCH.RNTI);
    
    % 36.212-5.3.3.3: Channel coding
    convblk = lteConvolutionalEncode(crcblk);

    % 36.212-5.3.3.4: Rate matching
    matchedblk = lteRateMatchConvolutional(convblk, ...
        72*2^enb.PDSCH.PDCCHFormat);

    % 36.212-6.8.2: Scrambling
    pdcchscrseq = ltePDCCHPRBS(enb,length(matchedblk));
    scrblk = xor(pdcchscrseq,matchedblk);
    
    % 36.212-6.8.3: Modulation
    modblk = lteSymbolModulate(scrblk,"QPSK");

    % Adds filler symbols if too few PDCCH symbols.
    pdcchblk = zeros(MTot/2,1);
    pdcchind = ltePDCCHSpace(enb,enb.PDSCH);
    pdcchblk((pdcchind(1,1)-1)/2+1:(pdcchind(1,1)-1)/2+length(modblk)) ...
        = modblk;

    % 6.8.4: Layer mapping and precoding
    layermapblk = lteLayerMap(pdcchblk,enb.PDSCH.NLayers);
    precodeblk = lteDLPrecode(enb,enb.PDSCH,layermapblk);

    % Multiplexing from 36.211-6.8.2
    multiplexblk = ltePDCCHInterleave(enb,precodeblk);
end

% Recovers CFI from PCFICH symbols.
function cfi = decodePCFICH(enb,pcfichsym)
    deprecodeblk = lteDLDeprecode(enb,enb.PDSCH,pcfichsym);
    layerdemapblk = lteLayerDemap(deprecodeblk,1);

    % demodblk are soft bits so unscrambling has to be done accordingly.
    demodblk = lteSymbolDemodulate(layerdemapblk{1},"QPSK");
    pdschscrseq = ltePCFICHPRBS(enb,length(demodblk));
    newsign = xor(pdschscrseq,sign(demodblk)==1);
    unscrblk = (newsign+-1*~newsign).*abs(demodblk);
    cfi = lteCFIDecode(unscrblk);
end

% Recovers DCI from PDCCH symbols.
function dci = decodePDCCH(enb,pdcchsym)
    % Computes length of DCI
    if enb.NDLRB <= 10
        P = 1;
    elseif enb.NDLRB <= 26
        P = 2;
    elseif enb.NDLRB <= 63
        P = 3;
    else
        P = 4;
    end

    dcibits = 14+ceil(enb.NDLRB/P);
    if enb.NDLRB <= 10
        dcibits = dcibits-1;
    end

    if dcibits == 15+ceil(log2(enb.NDLRB*(enb.NDLRB+1)/2))
        dcibits = dcibits+1;
    end

    if dcibits == 12 || dcibits == 14 || dcibits == 16 || dcibits == 20 ...
            || dcibits == 24 || dcibits == 26 || dcibits == 32 || ...
            dcibits == 40 || dcibits == 44 || dcibits == 56
        dcibits = dcibits+1;
    end
    recoverlen = 3*(dcibits+16);

    demultiplexblk = ltePDCCHDeinterleave(enb,pdcchsym);
    deprecodeblk = lteDLDeprecode(enb,enb.PDSCH,demultiplexblk);
    layerdemapblk = lteLayerDemap(deprecodeblk,1);

    % Removes filler symbols.
    pdcchind = ltePDCCHSpace(enb,enb.PDSCH);
    pdcchblk = layerdemapblk{1}((pdcchind(1,1)-1)/2+1: ...
        (pdcchind(1,1)-1)/2+72*2^enb.PDSCH.PDCCHFormat/2);

    % demodblk are soft bits so unscrambling has to be done accordingly.
    demodblk = lteSymbolDemodulate(pdcchblk,"QPSK");
    pdcchscrseq = ltePDCCHPRBS(enb,length(demodblk));
    newsign = xor(pdcchscrseq,sign(demodblk)==1);
    unscrblk = (newsign+-1*~newsign).*abs(demodblk);
    recoveredblk = lteRateRecoverConvolutional(unscrblk,recoverlen);
    deconvblk = lteConvolutionalDecode(recoveredblk);
    decrcblk = lteCRCDecode(deconvblk,"16",enb.PDSCH.RNTI);
    dci = lteDCI(enb,decrcblk);
end

% Recovers payload from PDSCH symbols
function decrcblk = decodePDSCH(enb,dci,pdschsym)
    % If MDS index is not recognized, stop decoding.
    if dci.ModCoding <= 9
        mod = "QPSK";
    elseif dci.ModCoding <= 16
        mod = "16QAM";
    else
        mod = "64QAM";
    end

    switch dci.ModCoding
        case 0
            translen = 152;
        case 4
            translen = 408;
        case 5
            translen = 4392;
        case 14
            translen = 12960;
        case 15
            translen = 14112;
        case 21
            translen = 6456;
        case 23
            translen = 12576;
        case 25
            if enb.NDLRB == 15
                translen = 8504;
            elseif enb.NDLRB == 25
                translen = 14112;
            else
                translen = 28336;
            end
        case 26
            if enb.NDLRB == 50
                translen = 30576;
            else
                translen = 45352;
            end
        case 27
            translen = 46888;
        otherwise
            decrcblk = -1;
            return;
    end

    enb.PDSCH.Modulation = mod;

    deprecodeblk = lteDLDeprecode(enb,enb.PDSCH,pdschsym);
    layerdemapblk = lteLayerDemap(deprecodeblk,1);

    % demodblk are soft bits so unscrambling has to be done accordingly.
    demodblk = lteSymbolDemodulate(layerdemapblk{1},mod);
    pdschscrseq = ltePDSCHPRBS(enb,enb.PDSCH.RNTI,0,length(demodblk));
    newsign = xor(pdschscrseq,sign(demodblk)==1);
    unscrblk = (newsign+-1*~newsign).*abs(demodblk);
    recoveredblk = lteRateRecoverTurbo(unscrblk,translen,dci.RV,enb.PDSCH);
    deturboblk = lteTurboDecode(recoveredblk);
    desegblk = lteCodeBlockDesegment(deturboblk,translen+24);
    decrcblk = lteCRCDecode(desegblk,"24A");
end