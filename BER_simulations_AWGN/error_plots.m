clear;clc;
load QPSK_ber_Eb_no.mat
load QPSK_ber_Es_no.mat
load QPSK_ser_Eb_no.mat
load QPSK_ser_Es_no.mat
load 4pam_gray.mat
load 4pam_uniform.mat
load 8psk_gray.mat
load 8psk_uniform.mat
load 16qam_gray.mat
load 16qam_uniform.mat
load bpsk.mat
load bfsk.mat
snr_db=0:2:20;
semilogy(snr_db, qpsk_ser_SNR_exp, 'r-x');hold on;grid on;
semilogy(snr_db, qpsk_ser_SNR_theory, 'rs');
semilogy(snr_db, qpsk_ser_Eb_No_exp, 'r-o');
semilogy(snr_db, qpsk_ser_Eb_No_theory, 'rd');
semilogy(snr_db, qpsk_ber_SNR_exp, 'b-x');
semilogy(snr_db, qpsk_ber_SNR_theory, 'bs');
semilogy(snr_db, qpsk_ber_Eb_No_exp, 'b-o');
semilogy(snr_db, qpsk_ber_Eb_No_theory, 'bd');
legend('SER E_s/N_o','SER E_s/N_o Theory','SER E_b/N_o','SER E_b/N_o Theory',...
    'BER E_s/N_o','BER E_s/N_o Theory','BER E_b/N_o','BER E_b/N_o Theory','Location','Southwest');
title('QPSK Error Analysis');
axis square;
set(gca,'FontSize',14);
xlabel('SNR (dB)');
ylabel('Error Rate/Probability');
ylim([1e-10 1]);
figure,
semilogy(snr_db, bpsk_ber_SNR_exp, 'r-x');hold on;grid on;
semilogy(snr_db, bpsk_ber_SNR_theory, 'rs');
legend('BER E_b/N_o','BER E_b/N_o Theory','Location','Southwest');
title('BPSK Error Analysis');
axis square;
set(gca,'FontSize',14);
xlabel('SNR (dB)');
ylabel('Error Rate/Probability');
ylim([1e-10 1]);
figure,
semilogy(snr_db, epsk_gray_ber_SNR_exp, 'b-x');hold on;grid on;
semilogy(snr_db, epsk_gray_ber_SNR_theory, 'bd');
semilogy(snr_db, epsk_uniform_ber_SNR_exp, 'r-x');hold on;grid on;
semilogy(snr_db, epsk_uniform_ber_SNR_theory, 'rs');
legend('BER Gray E_s/N_o','BER Gray E_s/N_o Theory','BER Uniform E_s/N_o','BER Uniform E_s/N_o Theory','Location','Southwest');
title('8-PSK Error Analysis');
axis square;
set(gca,'FontSize',14);
xlabel('SNR (dB)');
ylabel('Error Rate/Probability');
ylim([1e-10 1]);
figure,
semilogy(snr_db, pam_gray_ber_SNR_exp, 'b-x');hold on;grid on;
semilogy(snr_db, pam_gray_ber_SNR_theory, 'bd');
semilogy(snr_db, pam_uniform_ber_SNR_exp, 'r-x');hold on;grid on;
semilogy(snr_db, pam_uniform_ber_SNR_theory, 'rs');
legend('BER Gray E_s/N_o','BER Gray E_s/N_o Theory','BER Uniform E_s/N_o','BER Uniform E_s/N_o Theory','Location','Southwest');
title('4-PAM Error Analysis');
axis square;
set(gca,'FontSize',14);
xlabel('SNR (dB)');
ylabel('Error Rate/Probability');
ylim([1e-10 1]);
figure,
semilogy(snr_db, qam_gray_ber_SNR_exp, 'b-x');hold on;grid on;
semilogy(snr_db, qam_gray_ber_SNR_theory, 'bd');
semilogy(snr_db, qam_uniform_ber_SNR_exp, 'r-x');hold on;grid on;
semilogy(snr_db, qam_uniform_ber_SNR_theory, 'rs');
legend('BER Gray E_s/N_o','BER Gray E_s/N_o Theory','BER Uniform E_s/N_o','BER Uniform E_s/N_o Theory','Location','Southwest');
title('16-QAM Error Analysis');
axis square;
set(gca,'FontSize',14);
xlabel('SNR (dB)');
ylabel('Error Rate/Probability');
ylim([1e-10 1]);
figure,
semilogy(snr_db, bfsk_ber_SNR_exp, 'r-x');hold on;grid on;
semilogy(snr_db, bfsk_ber_SNR_theory, 'rs');
legend('BER E_b/N_o','BER E_b/N_o Theory','Location','Southwest');
title('BFSK Error Analysis');
axis square;
set(gca,'FontSize',14);
xlabel('SNR (dB)');
ylabel('Error Rate/Probability');
ylim([1e-10 1]);
figure,
semilogy(snr_db, bpsk_ber_SNR_exp, 'r-x');hold on;grid on;
semilogy(snr_db, bfsk_ber_SNR_exp, 'b-o');
semilogy(snr_db, epsk_gray_ber_SNR_exp, '-s');hold on;grid on;
semilogy(snr_db, pam_gray_ber_SNR_exp, '-*');
semilogy(snr_db, qam_gray_ber_SNR_exp, '-+');hold on;grid on;
semilogy(snr_db, qpsk_ber_SNR_exp, '-x');
legend('BPSK','BFSK','8-PSK','4-PAM','16-QAM','QPSK','Location','Southwest');
title('Comparison');
axis square;
set(gca,'FontSize',14);
xlabel('E_s/N_o(dB)');
ylabel('Error Rate/Probability');
ylim([1e-10 1]);