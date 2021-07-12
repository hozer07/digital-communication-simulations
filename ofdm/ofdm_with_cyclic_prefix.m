clear;clc;
qpsk_bits = [0 0;0 1;1 1;1 0];
qpsk_syms = exp((0:3).*1i.*pi./2+1i*pi/4);
N = 256;
num_blocks = 100;
channels{1} = [0.77 + 0.49*1i, 0.34 + 0.23*1i];
channels{2} = [0.20 + 0.74*1i, 0.22 + 0.36*1i, 0.18 - 0.29*1i, 0.08 - 0.27*1i, 0.07 - 0.20*1i];
SNR_range = 0:20;
ofdm_prefix = zeros(2,length(SNR_range));
err_found_lim = 100;
sim_lim = 1000;
for ch_index = [1,2]
    channel = channels{ch_index};
    channel_length = length(channel);
    for i=1:length(SNR_range)
        err_found_count = 0;
        sim_count = 0;
        No = 10^(-SNR_range(i)/10);
        while(err_found_count<err_found_lim && sim_count<sim_lim)
            sim_count = sim_count + 1;            
            data_syms = randi([1 4],1,N*num_blocks);            
            data_bits = qpsk_bits(data_syms,:);
            data_syms = qpsk_syms(data_syms);
            data_with_prefix = zeros(1,(N+channel_length-1)*num_blocks);
            for j=1:num_blocks
                data_syms((j-1)*N+1:j*N) = sqrt(N).*ifft(data_syms((j-1)*N+1:j*N));
                data_with_prefix((j-1)*(N+channel_length-1)+1:(j-1)*(N+channel_length-1)+channel_length-1) = data_syms(j*N-channel_length+2:j*N);
                data_with_prefix((j-1)*(N+channel_length-1)+channel_length:j*(N+channel_length-1))=data_syms((j-1)*N+1:j*N);                 
            end          
            received_syms = conv(data_with_prefix,channel);
            received_syms = received_syms + sqrt(No/2).*(randn(1,length(received_syms))+1i.*randn(1,length(received_syms)));
            for j=1:channel_length-1
                received_syms(1:(N+channel_length-j):end) = [];
            end
            received_syms2 = zeros(1,N*num_blocks);
            for j=1:num_blocks
                received_syms2((N*(j-1)+1):N*j) = 1/sqrt(N).*fft(received_syms((N*(j-1)+1):N*j))./fft(channel,N);
            end
            [~,received_syms2] = min(abs(repmat(received_syms2,4,1)-qpsk_syms.'));
            received_bits = qpsk_bits(received_syms2,:);
            err = sum(received_bits~=data_bits,'all');
            if(err>0),err_found_count=err_found_count+1;end
            ofdm_prefix(ch_index,i)=ofdm_prefix(ch_index,i)+err;
        end
        ofdm_prefix(ch_index,i) = ofdm_prefix(ch_index,i)/(sim_count*N*num_blocks*2);
    end
end
semilogy(SNR_range,ofdm_prefix(1,:),'*--','LineWidth',2,'MarkerSize',7);hold on;grid on;
semilogy(SNR_range,ofdm_prefix(2,:),'ro--','LineWidth',2,'MarkerSize',7);
% save part2a.mat ofdm_prefix
xlabel('E_s/N_0 (dB)');
ylabel('BER');
ylim([1e-5 1]);title('OFDM with cyclic prefix');
legend('Channel 1','Channel 2');