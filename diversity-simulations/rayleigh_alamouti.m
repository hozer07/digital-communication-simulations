clear;clc;
symbol_num = 1e5;
qpsk_syms = exp(1i.*(0:3)'.*pi/2+1i*pi/4);
SNR = 0:20;
sim_lim = 1000;
err_found_lim = 100;
fading_err_store = zeros(2,length(SNR));
for snr_index=1:length(SNR)
    noise_var = sqrt(10^(-SNR(snr_index)/10));
    %% 2 antennas MRC  
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
        received = received1.*channel1+received2.*channel2;
        [~,detected] = min(abs(repmat(received,1,4)-conj(qpsk_syms')),[],2);
        err = sum(detected~=qpsk_data,'all');
        if(err~=0),err_found_count=err_found_count+1;fading_err_store(1,snr_index)=fading_err_store(1,snr_index)+err;end        
    end
    fading_err_store(1,snr_index) = fading_err_store(1,snr_index)./sim_var./symbol_num;
    %% Alamouti
    sim_var = 0;
    err_found_count = 0;
    while(err_found_count<err_found_lim && sim_var<sim_lim)
        sim_var = sim_var+1;
        qpsk_data = randi([1 4],symbol_num,1);
        qpsk_symbols = qpsk_syms(qpsk_data);
        channel1 = sqrt(0.25).*(randn(1,symbol_num/2)+1i.*randn(1,symbol_num/2));
        channel2 = sqrt(0.25).*(randn(1,symbol_num/2)+1i.*randn(1,symbol_num/2));
        temp_data = reshape(qpsk_symbols,2,symbol_num/2);
        noise = noise_var/sqrt(2).*(randn(2,symbol_num/2)+1i.*randn(2,symbol_num/2));
        received_first_row = sum([channel1;channel2].*temp_data,1)+noise(1,:);
        received_second_row = sum([conj(channel2);-conj(channel1)].*temp_data,1)+noise(2,:);
        received = [received_first_row;received_second_row];
        received_first_row = sum([conj(channel1);channel2].*received,1);
        received_second_row = sum([conj(channel2);-channel1].*received,1);
        received = [received_first_row;received_second_row];
        received = reshape(received,1,symbol_num).';
        [~,detected] = min(abs(repmat(received,1,4)-qpsk_syms.'),[],2);
        err = sum(detected~=qpsk_data,'all');
        if(err~=0),err_found_count=err_found_count+1;fading_err_store(2,snr_index)=fading_err_store(2,snr_index)+err;end
    end
    fading_err_store(2,snr_index) = fading_err_store(2,snr_index)./sim_var./symbol_num;
end
semilogy(SNR,fading_err_store(1,:),'--*','MarkerSize',7);hold on;grid on;
semilogy(SNR,fading_err_store(2,:),'--+','MarkerSize',7);
ylim([1e-4 1]);
xlabel('SNR (dB)');ylabel('BER');
title('Rayleigh fading, 2 antenna MRC vs Alamouti');
legend('2 antenna MRC','Alamouti');