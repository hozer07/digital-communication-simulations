clear;clc;
symbol_num = 1e5;
qpsk_syms = exp(1i.*(0:3)'.*pi/2+1i*pi/4);
SNR = 0:20;
sim_lim = 1000;
err_found_lim = 100;
fading_err_store = zeros(4,length(SNR));
for snr_index=1:length(SNR)
    noise_var = sqrt(10^(-SNR(snr_index)/10));
    %% 1 antenna
    sim_var = 0;
    err_found_count = 0;
    while(err_found_count<err_found_lim && sim_var<sim_lim)
        sim_var = sim_var+1;
        qpsk_data = randi([1 4],symbol_num,1);
        qpsk_symbols = qpsk_syms(qpsk_data);
        channel = abs(sqrt(0.5).*(randn(symbol_num,1)+1i.*randn(symbol_num,1)));
        received = channel.*qpsk_symbols+noise_var/sqrt(2).*(randn(symbol_num,1)+1i.*randn(symbol_num,1));
        [~,detected] = min(abs(repmat(received./channel,1,4)-conj(qpsk_syms')),[],2);
        err = sum(detected~=qpsk_data);
        if(err~=0),err_found_count=err_found_count+1;fading_err_store(1,snr_index)=fading_err_store(1,snr_index)+err;end        
    end
    fading_err_store(1,snr_index) = fading_err_store(1,snr_index)./sim_var./symbol_num;
    %% 2 antennas    
    sim_var = 0;
    err_found_count = 0;
    while(err_found_count<err_found_lim && sim_var<sim_lim)
        sim_var = sim_var+1;
        qpsk_data = randi([1 4],symbol_num,1);
        qpsk_symbols = qpsk_syms(qpsk_data);
        channel1 = abs(sqrt(0.5).*(randn(symbol_num,1)+1i.*randn(symbol_num,1)));
        channel2 = abs(sqrt(0.5).*(randn(symbol_num,1)+1i.*randn(symbol_num,1)));
        received1 = (channel1).*qpsk_symbols+noise_var/sqrt(2).*(randn(symbol_num,1)+1i.*randn(symbol_num,1));
        received2 = (channel2).*qpsk_symbols+noise_var/sqrt(2).*(randn(symbol_num,1)+1i.*randn(symbol_num,1));
        received = channel1.*received1+channel2.*received2;
        [~,detected] = min(abs(repmat(received,1,4)-conj(qpsk_syms')),[],2);
        err = sum(detected~=qpsk_data,'all');
        if(err~=0),err_found_count=err_found_count+1;fading_err_store(2,snr_index)=fading_err_store(2,snr_index)+err;end        
    end
    fading_err_store(2,snr_index) = fading_err_store(2,snr_index)./sim_var./symbol_num;
    %% 3 antennas
    sim_var = 0;
    err_found_count = 0;
    while(err_found_count<err_found_lim && sim_var<sim_lim)
        sim_var = sim_var+1;
        qpsk_data = randi([1 4],symbol_num,1);
        qpsk_symbols = qpsk_syms(qpsk_data);
        channel1 = abs(sqrt(0.5).*(randn(symbol_num,1)+1i.*randn(symbol_num,1)));
        channel2 = abs(sqrt(0.5).*(randn(symbol_num,1)+1i.*randn(symbol_num,1)));
        channel3 = abs(sqrt(0.5).*(randn(symbol_num,1)+1i.*randn(symbol_num,1)));
        received1 = (channel1).*qpsk_symbols+noise_var/sqrt(2).*(randn(symbol_num,1)+1i.*randn(symbol_num,1));
        received2 = (channel2).*qpsk_symbols+noise_var/sqrt(2).*(randn(symbol_num,1)+1i.*randn(symbol_num,1));
        received3 = (channel3).*qpsk_symbols+noise_var/sqrt(2).*(randn(symbol_num,1)+1i.*randn(symbol_num,1));
        received = channel1.*received1+channel2.*received2+channel3.*received3;
        [~,detected] = min(abs(repmat(received,1,4)-conj(qpsk_syms')),[],2);
        err = sum(detected~=qpsk_data,'all');
        if(err~=0),err_found_count=err_found_count+1;fading_err_store(3,snr_index)=fading_err_store(3,snr_index)+err;end        
    end
    fading_err_store(3,snr_index) = fading_err_store(3,snr_index)./sim_var./symbol_num;
    %% 4 antennas
    sim_var = 0;
    err_found_count = 0;
    while(err_found_count<err_found_lim && sim_var<sim_lim)
        sim_var = sim_var+1;
        qpsk_data = randi([1 4],symbol_num,1);
        qpsk_symbols = qpsk_syms(qpsk_data);
        channel1 = abs(sqrt(0.5).*(randn(symbol_num,1)+1i.*randn(symbol_num,1)));
        channel2 = abs(sqrt(0.5).*(randn(symbol_num,1)+1i.*randn(symbol_num,1)));
        channel3 = abs(sqrt(0.5).*(randn(symbol_num,1)+1i.*randn(symbol_num,1)));
        channel4 = abs(sqrt(0.5).*(randn(symbol_num,1)+1i.*randn(symbol_num,1)));
        received1 = (channel1).*qpsk_symbols+noise_var/sqrt(2).*(randn(symbol_num,1)+1i.*randn(symbol_num,1));
        received2 = (channel2).*qpsk_symbols+noise_var/sqrt(2).*(randn(symbol_num,1)+1i.*randn(symbol_num,1));
        received3 = (channel3).*qpsk_symbols+noise_var/sqrt(2).*(randn(symbol_num,1)+1i.*randn(symbol_num,1));
        received4 = (channel4).*qpsk_symbols+noise_var/sqrt(2).*(randn(symbol_num,1)+1i.*randn(symbol_num,1));
        received = channel1.*received1+channel2.*received2+channel3.*received3+channel4.*received4;
        [~,detected] = min(abs(repmat(received,1,4)-conj(qpsk_syms')),[],2);
        err = sum(detected~=qpsk_data,'all');
        if(err~=0),err_found_count=err_found_count+1;fading_err_store(4,snr_index)=fading_err_store(4,snr_index)+err;end        
    end
    fading_err_store(4,snr_index) = fading_err_store(4,snr_index)./sim_var./symbol_num;
end
semilogy(SNR,fading_err_store(1,:),'--*','MarkerSize',7);hold on;grid on;
semilogy(SNR,fading_err_store(2,:),'--+','MarkerSize',7);
semilogy(SNR,fading_err_store(3,:),'--o','MarkerSize',7);
semilogy(SNR,fading_err_store(4,:),'--^','MarkerSize',7);ylim([1e-5 1]);
xlabel('SNR (dB)');ylabel('BER');
title('Rayleigh fading, QPSK, MRC');
legend('1 antenna','2 antennas','3 antennas','4 antennas');