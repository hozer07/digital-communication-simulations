clear;clc;
symbol_num = 1e5;
qpsk_syms = exp(1i.*(0:3)'.*pi/2+1i*pi/4);
SNR = 0:20;
sim_lim = 1000;
err_found_lim = 100;
alamouti_store = zeros(2,length(SNR));
spatial_mux = zeros(1,length(SNR));
for snr_index=1:length(SNR)
    noise_var = sqrt(10^(-SNR(snr_index)/10));
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
        if(err~=0),err_found_count=err_found_count+1;alamouti_store(1,snr_index)=alamouti_store(1,snr_index)+err;end
    end
    alamouti_store(1,snr_index) = alamouti_store(1,snr_index)./sim_var./symbol_num;
    %% Alamouti 2x2
    sim_var = 0;
    err_found_count = 0;
    while(err_found_count<err_found_lim && sim_var<sim_lim)
        sim_var = sim_var+1;
        qpsk_data = randi([1 4],symbol_num,1);
        qpsk_symbols = qpsk_syms(qpsk_data);
        channel_r1_t1 = sqrt(0.25).*(randn(1,symbol_num/2)+1i.*randn(1,symbol_num/2));
        channel_r1_t2 = sqrt(0.25).*(randn(1,symbol_num/2)+1i.*randn(1,symbol_num/2));
        channel_r2_t1 = sqrt(0.25).*(randn(1,symbol_num/2)+1i.*randn(1,symbol_num/2));
        channel_r2_t2 = sqrt(0.25).*(randn(1,symbol_num/2)+1i.*randn(1,symbol_num/2));
        temp_data = reshape(qpsk_symbols,2,symbol_num/2);
        noise = noise_var/sqrt(2).*(randn(4,symbol_num/2)+1i.*randn(4,symbol_num/2));
        received_first_row = sum([channel_r1_t1;channel_r1_t2].*temp_data,1)+noise(1,:);
        received_second_row = sum([-conj(channel_r1_t2);conj(channel_r1_t1)].*temp_data,1)+noise(2,:);
        received_third_row = sum([channel_r2_t1;channel_r2_t2].*temp_data,1)+noise(3,:);
        received_fourth_row = sum([-conj(channel_r2_t2);conj(channel_r2_t1)].*temp_data,1)+noise(4,:);
        received = [received_first_row;received_second_row;received_third_row;received_fourth_row];
        received_first_row = sum([conj(channel_r1_t1);-channel_r1_t2;conj(channel_r2_t1);-channel_r2_t2].*received,1);
        received_second_row = sum([conj(channel_r1_t2);channel_r1_t1;conj(channel_r2_t2);channel_r2_t1].*received,1);
        received = reshape([received_first_row;received_second_row],1,symbol_num);
        [~,detected] = min(abs(repmat(received.',1,4)-qpsk_syms.'),[],2);
        err = sum(detected~=qpsk_data,'all');
        if(err~=0),err_found_count=err_found_count+1;alamouti_store(2,snr_index)=alamouti_store(2,snr_index)+err;end
    end
    alamouti_store(2,snr_index) = alamouti_store(2,snr_index)./sim_var./symbol_num;
    %% Alamouti 2x2 Spatial Mux
%     sim_var = 0;
%     err_found_count = 0;
%     while(err_found_count<err_found_lim && sim_var<sim_lim)
%         sim_var = sim_var+1;
%         qpsk_data = randi([1 4],symbol_num,1);
%         qpsk_symbols = qpsk_syms(qpsk_data);
%         channel_r1_t1 = sqrt(0.25).*(randn(1,symbol_num)+1i.*randn(1,symbol_num));
%         channel_r1_t2 = sqrt(0.25).*(randn(1,symbol_num)+1i.*randn(1,symbol_num));
%         channel_r2_t1 = sqrt(0.25).*(randn(1,symbol_num)+1i.*randn(1,symbol_num));
%         channel_r2_t2 = sqrt(0.25).*(randn(1,symbol_num)+1i.*randn(1,symbol_num));      
%         temp_data = reshape(qpsk_symbols,2,symbol_num/2);
%         noise = noise_var/sqrt(2).*(randn(2,symbol_num/2)+1i.*randn(2,symbol_num/2));
%         received = zeros(symbol_num,1);
%         for i=1:symbol_num/2
%             [U,S,V] = svd([channel_r1_t1(i) channel_r1_t2(i);channel_r2_t1(i) channel_r2_t2(i)]);
%             precoded = V*temp_data(:,i);
%             received(2*i-1:2*i) = [channel_r1_t1(i) channel_r1_t2(i);channel_r2_t1(i) channel_r2_t2(i)]*precoded+noise(:,i);
%             received(2*i-1:2*i) = U'*received(2*i-1:2*i);
%         end
%         [~,detected] = min(abs(repmat(received,1,4)-qpsk_syms.'),[],2);
%         err = sum(detected~=qpsk_data,'all');
%         if(err~=0),err_found_count=err_found_count+1;spatial_mux(1,snr_index)=spatial_mux(1,snr_index)+err;end
%     end
%     spatial_mux(1,snr_index) = spatial_mux(1,snr_index)./sim_var./symbol_num;
end
load spatial_mux.mat
semilogy(SNR,alamouti_store(1,:),'--*','MarkerSize',7);hold on;grid on;
semilogy(SNR,alamouti_store(2,:),'--+','MarkerSize',7);
semilogy(SNR,spatial_mux,'--+','MarkerSize',7);
ylim([1e-4 1]);
xlabel('SNR (dB)');ylabel('BER');
legend('Alamouti 2x1','Alamouti 2x2','Alamouti Spatial Multiplexing');