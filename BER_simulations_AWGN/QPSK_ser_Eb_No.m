clear all
warning off
%%%%% SYMBOL CONSTELLATION %%%%%%%%%
tb4=[0 0;0 1; 1 1; 1 0];
tb8=[zeros(4,1) tb4; ones(4,1) tb4];
tb16=[zeros(8,1) tb8; ones(8,1) tb8];
nBitPerSym=2;
M=2^nBitPerSym;
symbolBook=exp(1j*2*pi*[0:(M-1)]/M);
bitBook=tb4;
%%%%%%%%%%%%%% ERROR MATRICES %%%%%%%%%%%%%%%%%%%%%%%
nBitsPerFrame=3000;
max_nFrame=2000;
fErrLim=100;
snr_db=0:2:20;
nSymPerFrame=nBitsPerFrame/nBitPerSym;
errs=zeros(length(snr_db), 1);
nFrames=zeros(length(snr_db), 1);
fErrs=errs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nEN = 1:length(snr_db) % SNR POINTS
    snr_p=snr_db(nEN);
    en = 10^(snr_p/10); % convert SNR from unit db to normal numbers
    sigma = 1/sqrt(en); % standard deviation of AWGN noise
    nframe = 0;
    while (nframe<max_nFrame) && (fErrs(nEN)<fErrLim)
        err_count=0;
        nframe = nframe + 1;
        info_bits=round(rand(1,nBitsPerFrame));
        info_matrix=reshape(info_bits, nBitPerSym, nSymPerFrame);
        sym_vec=ones(1,nSymPerFrame);
        for v=1:nSymPerFrame
            for i=1:M
                if(info_matrix(:,v)' == bitBook(i,:))
                    sym_vec(v)=i;
                    break;
                end
            end
        end
        sym_seq=symbolBook(sym_vec);
        %%%%%%%%%%%%%CHANNEL %%%%%%%%%%%%%%%%%%%%%
        noise=1/2*[randn(1, nSymPerFrame) + 1j*randn(1,nSymPerFrame)];
        det_seq=zeros(1,nBitsPerFrame);
        rec_sig=sym_seq+sigma*noise;
        %%%%DETECTOR %%%%%%%%%%%%
        symbolBook_mat=repmat(transpose(symbolBook),1,nSymPerFrame);
        rec_sig_mat=repmat(rec_sig,length(symbolBook),1);
        distance_mat=abs(symbolBook_mat-rec_sig_mat);
        [~, det_sym_ind]=min(distance_mat,[],1);
        detected_bits=[bitBook(det_sym_ind, :)]';
        err = sum(round(sum(abs(info_matrix-detected_bits))./2));
        errs(nEN)=errs(nEN)+err;
        err_count=err_count+err;
        if err_count~=0
            fErrs(nEN)=fErrs(nEN)+1;
        end
    end % End of while loop
    nFrames(nEN)=nframe;
    sim_res=[errs nFrames]
end %end for (SNR points)
SNR = 10.^(snr_db./10);
Eb_No = SNR.*nBitPerSym;
qpsk_ser_Eb_No_exp = errs./nFrames/nSymPerFrame;
qpsk_ser_Eb_No_theory = 2.*qfunc(sqrt(2.*SNR)).*(1-0.5.*qfunc(sqrt(2.*SNR)));
save QPSK_ser_Eb_no.mat qpsk_ser_Eb_No_exp qpsk_ser_Eb_No_theory