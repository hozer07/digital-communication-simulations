% Frequency domain equalization with 1,3,5,10
% bits long cyclic prefixes for BPSK are simulated in this simulation.
clear;
clc;
snr_start = 0;
snr_increment = 2;
Eb_No_dB = snr_start:snr_increment:20;
cyclic_prefix_lengths = [1 3 5 10];  % Cyclic prefix lengths in the simulation
FDEsimBer = zeros(length(Eb_No_dB),length(cyclic_prefix_lengths)); % Placeholder for results
sim_count = zeros(size(FDEsimBer)); % Placeholder for simulation counts
block_size = 1000; % Size of an information block
channel = [0.74 -0.514 0.37 0.216 0.062]; % Impulse response of the channel
err_lim = 500; % Minimum number of errors for an SNR level and prefix length
err_check = false; % Flag to check if "err_lim" is satisfied
for snr=1:length(Eb_No_dB)
    for cyclic=1:length(cyclic_prefix_lengths)
        noise_variance = 1/(2*10^(Eb_No_dB(snr)/10));
        sim_num = 0;
        while(sim_num < 1000 || not(err_check))
            sim_num = sim_num+1;
            data = 2.*randi([0 1],1,block_size+cyclic_prefix_lengths(cyclic))-1; % Uniformly generated BPSK data
            data(1:cyclic_prefix_lengths(cyclic))=data(end-cyclic_prefix_lengths(cyclic)+1:end); % Prefix is added
            noise = sqrt(noise_variance).*randn(1,length(data)+length(channel)-1);
            y = conv(data,channel)+noise; % Linear convolution of channel and data      
            y_out = y(cyclic_prefix_lengths(cyclic)+1:cyclic_prefix_lengths(cyclic)+block_size); % Prefix is discarded
            equalized = ifft(fft(y_out)./fft([channel zeros(1,length(y_out)-length(channel))])); % Equalization assuming circular convolution
            equalized(equalized>0)=1;equalized(equalized<0)=-1; % ML detection
            FDEsimBer(snr,cyclic) = FDEsimBer(snr,cyclic)+sum(equalized~=data(cyclic_prefix_lengths(cyclic)+1:end));
            if(FDEsimBer(snr,cyclic)>=err_lim),err_check = true;end
        end
        sim_count(snr,cyclic) = sim_num;
        err_check = false;
    end
end
figure
semilogy(Eb_No_dB,FDEsimBer(:,1)./block_size./sim_count(:,1),'r-','Linewidth',1);
grid on;hold on;
semilogy(Eb_No_dB,FDEsimBer(:,2)./block_size./sim_count(:,2),'g-','Linewidth',1);
semilogy(Eb_No_dB,FDEsimBer(:,3)./block_size./sim_count(:,3),'b-','Linewidth',1);
semilogy(Eb_No_dB,FDEsimBer(:,4)./block_size./sim_count(:,4),'m-','Linewidth',1);
legend('1 prefix','3 prefix','5 prefix','10 prefix')
xlabel('Eb/No dB');
ylabel('Bit Error Rate');
axis square;