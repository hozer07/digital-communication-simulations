clear;clc;
symbol_num = 1e5;
bpsk_syms = [-1 1];
qpsk_syms = exp(1i.*(0:3)'.*pi/2+1i*pi/4);
qpsk_bits = [0 0;0 1;1 1;1 0];
eight_psk_syms = exp(1i.*(0:7)'.*pi/4);
eight_psk_bits = [zeros(4,1) qpsk_bits; ones(4,1) flipud(qpsk_bits)];
sixteen_qam_syms = qammod(0:15,16)./sqrt(10);
sixteen_qam_bits = [zeros(8,1) eight_psk_bits; ones(8,1) flipud(eight_psk_bits)];
SNR = 0:20;
sim_lim = 1000;
err_found_lim = 100;
fading_err_store = zeros(4,length(SNR));
awgn_err_store = zeros(4,length(SNR));
pd = makedist('Nakagami','mu',2,'omega',1);
for snr_index=1:length(SNR)
    noise_var = sqrt(10^(-SNR(snr_index)/10));
    %% BPSK simulation
    sim_var = 0;
    err_found_count = 0;
    while(err_found_count<err_found_lim && sim_var<sim_lim)
        sim_var = sim_var+1;
        bpsk_data = randi([1 2],1,symbol_num);
        bpsk_symbols = bpsk_syms(bpsk_data);
        channel = random(pd,1,symbol_num);
        received = channel.*bpsk_symbols+noise_var/sqrt(2)*(randn(1,symbol_num));
        detected = (received./channel>0)+1;
        err = sum(detected~=bpsk_data);
        if(err~=0),err_found_count=err_found_count+1;fading_err_store(1,snr_index)=fading_err_store(1,snr_index)+err;end        
    end
    fading_err_store(1,snr_index) = fading_err_store(1,snr_index)./sim_var./symbol_num;
    %% QPSK simulation
    sim_var = 0;
    err_found_count = 0;
    while(err_found_count<err_found_lim && sim_var<sim_lim)
        sim_var = sim_var+1;
        qpsk_data = randi([1 4],symbol_num,1);
        qpsk_bit_data = qpsk_bits(qpsk_data,:);
        qpsk_symbols = qpsk_syms(qpsk_data);
        channel = random(pd,symbol_num,1);
        received = channel.*qpsk_symbols+noise_var/sqrt(2).*(randn(symbol_num,1)+1i.*randn(symbol_num,1));
        [~,detected] = min(abs(repmat(received./channel,1,4)-conj(qpsk_syms')),[],2);
        detected = qpsk_bits(detected,:);
        err = sum(detected~=qpsk_bit_data,'all');
        if(err~=0),err_found_count=err_found_count+1;fading_err_store(2,snr_index)=fading_err_store(2,snr_index)+err;end        
    end
    fading_err_store(2,snr_index) = fading_err_store(2,snr_index)./sim_var./symbol_num./2;
    %% 8-PSK simulation
    sim_var = 0;
    err_found_count = 0;
    while(err_found_count<err_found_lim && sim_var<sim_lim)
        sim_var = sim_var+1;
        eight_psk_data = randi([1 8],symbol_num,1);
        eight_bit_data = eight_psk_bits(eight_psk_data,:);
        eight_psk_symbols = eight_psk_syms(eight_psk_data);
        channel = random(pd,symbol_num,1);
        received = channel.*eight_psk_symbols+noise_var/sqrt(2).*(randn(symbol_num,1)+1i.*randn(symbol_num,1));
        [~,detected] = min(abs(repmat(received./channel,1,8)-conj(eight_psk_syms')),[],2);
        detected = eight_psk_bits(detected,:);
        err = sum(detected~=eight_bit_data,'all');
        if(err~=0),err_found_count=err_found_count+1;fading_err_store(3,snr_index)=fading_err_store(3,snr_index)+err;end        
    end
    fading_err_store(3,snr_index) = fading_err_store(3,snr_index)./sim_var./symbol_num./3;
    %% 16-QAM simulation
    sim_var = 0;
    err_found_count = 0;
    while(err_found_count<err_found_lim && sim_var<sim_lim)
        sim_var = sim_var+1;
        sixteen_qam_data = randi([1 16],symbol_num,1);
        sixteen_qam_bit_data = sixteen_qam_bits(sixteen_qam_data,:);
        sixteen_qam_symbols = sixteen_qam_syms(sixteen_qam_data).';
        channel = random(pd,symbol_num,1);
        received = channel.*sixteen_qam_symbols+noise_var/sqrt(2).*(randn(symbol_num,1)+1i.*randn(symbol_num,1));
        [~,detected] = min(abs(repmat(received./channel,1,16)-sixteen_qam_syms),[],2);
        detected = sixteen_qam_bits(detected,:);
        err = sum(detected~=sixteen_qam_bit_data,'all');
        if(err~=0),err_found_count=err_found_count+1;fading_err_store(4,snr_index)=fading_err_store(4,snr_index)+err;end        
    end
    fading_err_store(4,snr_index) = fading_err_store(4,snr_index)./sim_var./symbol_num./4;
end
step_size = pi/20000;
interval = 0:step_size:pi/2;
mgf = @(beta,theta,SNR_dec) (1+beta./sin(theta).^2.*SNR_dec./4).^(-2);
for i=1:length(SNR)
    SNR_dec = 10^(SNR(i)/10);
    theoretical_bpsk(i) = 1/pi.*sum(mgf(2,interval,SNR_dec)).*step_size;
    theoretical_qpsk(i) = 1/pi.*sum(mgf(1,interval,SNR_dec)).*step_size;
    theoretical_8psk(i) = 2/3/pi.*sum(mgf(2*sin(pi/8)^2,interval,SNR_dec)).*step_size;
    theoretical_16qam(i) = 1/pi.*sum(mgf(0.2,interval,SNR_dec)).*step_size;
end
semilogy(SNR,fading_err_store(1,:),'-*','MarkerSize',8);hold on;
semilogy(SNR,fading_err_store(2,:),'--*','LineWidth',2);
semilogy(SNR,fading_err_store(3,:),'-*','LineWidth',2);
semilogy(SNR,fading_err_store(4,:),'-*','LineWidth',2);
semilogy(SNR,theoretical_bpsk,'ro','MarkerSize',8);
semilogy(SNR,theoretical_qpsk,'o','MarkerSize',8);
semilogy(SNR,theoretical_8psk,'^','MarkerSize',8);
semilogy(SNR,theoretical_16qam,'d','MarkerSize',8);
legend('BPSK simulation','QPSK simulation','8-PSK simulation','16-QAM simulation',...
    'Theoretical BPSK','Theoretical QPSK','Theoretical 8-PSK','Theoretical 16-QAM');
axis square;
grid on;
xlabel ('E_s/N_0 (dB)');
ylabel('BER');
title('Nakagami-m fading m=2');