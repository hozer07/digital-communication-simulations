% IMPORTANT
% The difference between theoretical symbol error probability and
% simulation symbol error rate for Rayleigh fading channel is due to
% the assumption that symbol errors are caused by only 1 bit error;
% however, 2 bit errors are also present in the simulation, though they are
% very rare. As you will see, BER and BEP plots exactly overlap,
% but SER and SEP plots have a narrow gap in between.
clear;clc;
err_min=100;
err_check = true;
sim_lim = 20;                               % Number of simulations to average over for a given SNR
qpsk_sym = exp(1i*pi.*(1:2:7)./4);          % Baseband QPSK symbols
qpsk_bits =[0,0;0,1;1,1;1,0];               % Gray coded bits for each symbol
SNR_increment = 5;                          % SNR increment for simulations
SNR_start = 0;                              % Minimum SNR
SNR_int = SNR_start:SNR_increment:25;       % SNR range
fading_symbol_errs_store = zeros(length(SNR_int),2);      % Placeholder matrix for symbol errors in fading case
awgn_ser_store = zeros(length(SNR_int),1);        % Placeholder matrix for symbol errors in awgn-only case
awgn_ber_store = zeros(length(SNR_int),1);        % Placeholder matrix for bit errors in awgn-only case
fading_ber_store = zeros(length(SNR_int),1);      % Placeholder matrix for bit errors in fading case
seq_len = 1e5;                              % Number of symbols for each simulation
for snr=SNR_int
    snr_index = snr/SNR_increment-SNR_start/SNR_increment+1;
    noise_stdev = sqrt(1/(2*2*10^(snr/10)));    % sqrt(No/2)
    sim_num = 0;
    while (err_check==true || sim_num <= sim_lim)
        % Rayleigh Fading channel
        sim_num = sim_num + 1;
        info_seq = randi([0,3],seq_len,1);      % Uniformly generated symbols
        info_bits = qpsk_bits(info_seq+1,:);    % Bits for generated symbols
        info_bits2 = info_bits';
        info_bits3 = info_bits2(:);
        channel = sqrt(0.5).*(randn(seq_len,1)+1i.*randn(seq_len,1));   % Rayleigh fading channel
        four_sym_seq1 = conj(qpsk_sym(info_seq+1)');
        four_sym_seq = four_sym_seq1.*abs(channel);                      % Symbols multiplied by channel gain
        four_sym_seq = four_sym_seq + noise_stdev.*(randn(seq_len,1)+1i.*randn(seq_len,1)); % Noise added
        four_sym_seq = repmat(four_sym_seq,[1,4]);
        fading_norms = abs(four_sym_seq-repmat(qpsk_sym,[seq_len,1]));
        [~,fading_detected] = min(fading_norms,[],2);           % Maximum Likelihood detection
        fading_detected = fading_detected-1;
        fading_symbol_err = length(find(fading_detected~=info_seq));
        fading_detected_bits = qpsk_bits(fading_detected+1,:);
        fading_detected_bits2 = fading_detected_bits';
        fading_detected_bits3 = fading_detected_bits2(:);
        fading_ber = length(find(fading_detected_bits3~=info_bits3));
        fading_ber_store(snr_index) = fading_ber_store(snr_index) + fading_ber;
        fading_symbol_errs_store(snr_index,1)=fading_symbol_errs_store(snr_index,1)+fading_symbol_err;
        % AWGN-only channel
        awgn_four_sym_seq = four_sym_seq1 + noise_stdev.*(randn(seq_len,1)+1i.*randn(seq_len,1));
        awgn_norms = abs(awgn_four_sym_seq-repmat(qpsk_sym,[seq_len,1]));
        [~,awgn_detected] = min(awgn_norms,[],2);
        awgn_detected = awgn_detected-1;
        awgn_ser = length(find(awgn_detected~=info_seq));
        awgn_detected_bits = qpsk_bits(awgn_detected+1,:);
        awgn_detected_bits2 = awgn_detected_bits';
        awgn_detected_bits3 = awgn_detected_bits2(:);
        awgn_ber = length(find(awgn_detected_bits3~=info_bits3));
        awgn_ber_store(snr_index) = awgn_ber_store(snr_index) + awgn_ber;
        awgn_ser_store(snr_index)=awgn_ser_store(snr_index,1)+awgn_ser;      
        if(fading_symbol_errs_store(snr_index,1)>100),err_check=false;end
    end
    fading_symbol_errs_store(snr_index,2)=sim_num;
end
err_store = fading_symbol_errs_store(:,1)./fading_symbol_errs_store(:,2)./seq_len;  % Averaging the error probabilities by
fading_ber_store = fading_ber_store./fading_symbol_errs_store(:,2)./seq_len./2;     % dividing total number of errors
awgn_ser_store = awgn_ser_store./seq_len./fading_symbol_errs_store(:,2);            % by simulation count
awgn_ber_store = awgn_ber_store./seq_len./fading_symbol_errs_store(:,2)./2;
EbNo = 10.^(SNR_int./10);
theoretical_symbol_error_qpsk_noise_only = 2.*qfunc(sqrt(2.*EbNo)).*(1-0.5.*qfunc(sqrt(2.*EbNo)));
theoretical_bit_error_qpsk_noise_only = theoretical_symbol_error_qpsk_noise_only./2;
theoretical_bit_error_qpsk_rayleigh = 0.5.*(1-sqrt(EbNo./(1+EbNo)));
theoretical_symbol_error_qpsk_rayleigh = 2.*theoretical_bit_error_qpsk_rayleigh;
semilogy(SNR_int,err_store,'-*');hold on;grid on;
semilogy(SNR_int,theoretical_symbol_error_qpsk_rayleigh,'-d');
semilogy(SNR_int,awgn_ser_store,'-x','LineWidth',1.7);
semilogy(SNR_int,theoretical_symbol_error_qpsk_noise_only,'-s');
legend('Rayleigh Simulation SER','Rayleigh Theoretical SER','AWGN-only Simulation','AWGN-only Theoretical');xlabel('E_b/N_o');
title('SER Performance of QPSK in Rayleigh Fading Channel')
set(gca,'FontSize',10);ylim([1e-6 1]);
figure,
semilogy(SNR_int,fading_ber_store,'-*');hold on;grid on;
semilogy(SNR_int,theoretical_bit_error_qpsk_rayleigh,'-d');
semilogy(SNR_int,awgn_ber_store,'-x','LineWidth',1.7);ylim([1e-6 1]);
semilogy(SNR_int,theoretical_bit_error_qpsk_noise_only,'-s');
legend('Rayleigh Simulation BER','Rayleigh Theoretical BER','AWGN-only Simulation','AWGN-only Theoretical');
xlabel('E_b/N_o');title('BER Performance Of QPSK in Rayleigh Fading Channel')
set(gca,'FontSize',10);ylim([1e-6 1]);
