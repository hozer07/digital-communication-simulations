% Three receiver diversity methods; selection diversity,
% maximum ratio combining and equal gain combining
% for 2 Rayleigh fading channels are
% compared in terms of error probabilities
% in this simulation.
clear;clc;
err_min=100;
err_check = true;
sim_lim = 20;
qpsk_sym = pskmod(0:3,4);                   % Baseband QPSK symbols
qpsk_bits =[0,0;0,1;1,1;1,0];               % Gray coded bits for each symbol
SNR_start = 0;
SNR_increment = 3;
SNR_int = SNR_start:SNR_increment:25;
errs_store = zeros(3,length(SNR_int),2);
ber_div = zeros(3,length(SNR_int),1);
seq_len = 1e5;                              % Number of symbols for each simulation
for method = 1:3                    
    for snr=SNR_int
        snr_index = snr/SNR_increment-SNR_start/SNR_increment+1;
        noise_stdev = sqrt(1/(2*2*10^(snr/10)));
        sim_num = 0;
        while (err_check==true || sim_num <= sim_lim)
            sim_num = sim_num + 1;
            info_seq = randi([0,3],seq_len,1);
            info_bits = qpsk_bits(info_seq+1,:);
            info_bits2 = info_bits';
            info_bits3 = info_bits2(:);
            channel1 = sqrt(0.25).*(randn(seq_len,1)+1i.*randn(seq_len,1)); % multiplication by sqrt(0.25) is to
            channel2 = sqrt(0.25).*(randn(seq_len,1)+1i.*randn(seq_len,1)); % keep overall variance 1
            four_sym_seq = conj(qpsk_sym(info_seq+1)');
            four_sym_seq1 = four_sym_seq.*abs(channel1);
            four_sym_seq2 = four_sym_seq.*abs(channel2);
            four_sym_seq1 = four_sym_seq1 + noise_stdev.*(randn(seq_len,1)+1i.*randn(seq_len,1));
            four_sym_seq2 = four_sym_seq2 + noise_stdev.*(randn(seq_len,1)+1i.*randn(seq_len,1));
            if (method==1) % Selection diversity
                [~,indices] = max(abs([four_sym_seq1 four_sym_seq2]),[],2);
                four_sym_seq = four_sym_seq1.*(indices==1)+four_sym_seq2.*(indices==2);
            elseif(method==2),four_sym_seq=(four_sym_seq1+four_sym_seq2); % Equal gain combining
            else,four_sym_seq = abs(channel1).*four_sym_seq1+abs(channel2).*four_sym_seq2; % Maximum ratio combining
            end
            four_sym_seq = repmat(four_sym_seq,[1,4]);
            norms = abs(four_sym_seq-repmat(qpsk_sym,[seq_len,1]));
            [~,detected] = min(norms,[],2); % Maximum likelihood detection
            detected = detected-1;
            err = length(find(detected~=info_seq));
            detected_bits = qpsk_bits(detected+1,:);
            detected_bits2 = detected_bits';
            detected_bits3 = detected_bits2(:);
            ber_temp = length(find(detected_bits3~=info_bits3));
            ber_div(method,snr_index) = ber_div(method,snr_index) + ber_temp;
            errs_store(method,snr_index,1)=errs_store(method,snr_index,1)+err;
            if(errs_store(method,snr_index,1)>100),err_check=false;end
        end
        errs_store(method,snr_index,2)=sim_num;
    end
 end
err_store_div = errs_store(:,:,1)./errs_store(:,:,2)./seq_len;
ber_div = ber_div./errs_store(:,:,2)./seq_len./2;
EbNo = 10.^(SNR_int./10);
load no_diversity.mat
semilogy(SNR_int,err_store,'-*');hold on;grid on;
semilogy(SNR_int,err_store_div(1,:),'-*');hold on;grid on;
semilogy(SNR_int,err_store_div(2,:),'-x','LineWidth',1.7);
semilogy(SNR_int,err_store_div(3,:),'-','LineWidth',1.7);legend('No diversity','Selection Diversity',...
    'Equal Gain Combining','Maximal Ratio Combining');xlabel('E_b/N_o');ylabel('SER');title('SER comparison of diversity methods');
set(gca,'FontSize',10);ylim([1e-5 1]);
figure,
semilogy(SNR_int,ber,'-*');hold on;grid on;
semilogy(SNR_int,ber_div(1,:),'-*');hold on;grid on;
semilogy(SNR_int,ber_div(2,:),'-x','LineWidth',1.7);
semilogy(SNR_int,ber_div(3,:),'-','LineWidth',1.7);legend('No diversity','Selection Diversity',...
    'Equal Gain Combining','Maximal Ratio Combining');xlabel('E_b/N_o');ylabel('BER');title('BER comparison of diversity methods');
set(gca,'FontSize',10);ylim([1e-5 1]);