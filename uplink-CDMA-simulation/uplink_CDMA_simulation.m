% This simulation demonstrates the BER performance
% of uplink CDMA for 2,4 and 8 user cases under Rayleigh fading channels.
clear;clc;
sim_lim = 20; % Number of simulations to average over for each SNR level
snr_start = -10; %dB
snr_increment = 5;
SNR_int = snr_start:snr_increment:20;
number_of_users = [2,4,8]; % Number of users talking at the same time
errs_store = zeros(length(number_of_users),length(SNR_int)); % Placeholder matrix for results
seq_len = 1e4; % number of bits for each simualtion
reg = 4; %LFSR register size, LFSR generator supports 4 or 5
user = 1; %determines whose despreading code will be used, should be in the range (1,number_of_users)
code_len = 2^reg-1;  % Length of the spreading code
[~,codes] = lfsr_generator(reg,max(number_of_users));
for K=1:length(number_of_users)
    for snr=SNR_int
        snr_index = snr/snr_increment-snr_start/snr_increment+1;
        noise_stdev = sqrt(1/(2*10^(snr/10)));
        sim_num = 0;
        while (sim_num <= sim_lim)
            sim_num = sim_num + 1;
            channels = abs(sqrt(0.5).*(randn(number_of_users(K),seq_len*code_len)+1i.*randn(number_of_users(K),seq_len*code_len))); % Rayleigh
%             channels = ones(number_of_users(K),seq_len*code_len); %Uncomment this line for non-fading channels
            noise = noise_stdev.*(randn(number_of_users(K),seq_len*code_len));
            info_bits = randi([0 1],[number_of_users(K) seq_len]);        
            % spreading by xor operation and BPSK symbol conversion.
            spreaded_seq = mod(repelem(info_bits,1,code_len)+repmat(codes(1:number_of_users(K),:),[1 seq_len]),2);
            spreaded_seq(spreaded_seq==1) = -1;spreaded_seq(spreaded_seq==0) = 1;
            received = sum(spreaded_seq.*channels+noise);
            % despreading in the same way as spreading and detection
            despreading_code = codes(user,:);despreading_code(despreading_code==1)= -1;despreading_code(despreading_code==0)= 1;
            despreaded_seq = received.*repmat(despreading_code,[1 seq_len]);
            detected = cumsum(despreaded_seq); detected = detected(code_len:code_len:end);detected = [detected(1) diff(detected)]./code_len;
            final_detected=detected;final_detected(detected<0)=1;final_detected(detected>0)=0;
            err = length(find(final_detected~=info_bits(user,:)));
            errs_store(K,snr_index)=errs_store(K,snr_index)+err;
        end
    end
end
% Add or remove plot lines if you want to examine 
% different number of users.
err_store = errs_store./sim_lim./seq_len;
EbNo = 10.^(SNR_int./10);
semilogy(SNR_int,err_store(1,:),'-*');hold on;grid on;
semilogy(SNR_int,err_store(2,:),'-+');hold on;grid on;
semilogy(SNR_int,err_store(3,:),'-^');hold on;grid on;
xlabel('E_b/N_o');ylabel('BER');title('BER for m=4 codes, flat fading channel');
legend('K=2','K=4','K=8');
set(gca,'FontSize',10);