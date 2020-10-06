clear;clc
N  = 10^5; % number of bits to calculate BER
SNR_start = 0;
SNR_increment = 2;
tap_options = [9 15]; % Lengths of the equalizers in the simulation
Eb_No_dB = SNR_start:SNR_increment:40; % Eb/No Range
%Impulse Responses of the channels in the simulation
channels{1} = [0.04 -0.06 0.07 -0.21 -0.5 0.72 0.36 0.21 0.06];
channels{2} = [0.41 0.815 0.41];
channels{3} = [0.227 0.460 0.688 0.460 0.227];
MMSEsimBer = zeros(length(Eb_No_dB),length(channels),length(tap_options)); % Placeholder for simulation results
for average = 1:100 % The number of simulations to average the results
    for snr=1:length(Eb_No_dB)
        infos = 2.*randi([0 1],1,N)-1; % Uniformly generated binary information sequence
        noise_variance = 0.5*(10^(-Eb_No_dB(snr)/10));
        noise = sqrt(noise_variance).*randn(1,N);
        for channel_index = 1:length(channels)
            channel_output = conv(infos,channels{channel_index},'same');
            y = channel_output+noise;
            system_length = length(channels{channel_index});
            for tap=1:length(tap_options)
                % c = argmin|conv(c,y)-infos|, satisfying c is
                % found below
                channel_matrix = toeplitz([channels{channel_index}(1);zeros(tap_options(tap)-1,1)],[channels{channel_index} zeros(1,tap_options(tap)-1)]);          
                b = zeros(system_length+tap_options(tap)-1,1);
                b(round((system_length+tap_options(tap)-1)/2)) = 1;
                c = (channel_matrix*channel_matrix'+noise_variance.*eye(size(channel_matrix,1)))\channel_matrix*b;           
                y_out_ZF = conv(y,c,'same');
                y_out_ZF = 2.*double(y_out_ZF>0)-1;
                MMSEsimBer(snr,channel_index,tap) = MMSEsimBer(snr,channel_index,tap)+size(find(y_out_ZF~=infos),2);%Store zero forcing results         
            end
        end
    end
end
MMSEsimBer = MMSEsimBer./N./average; % simulated ber
figure
semilogy(Eb_No_dB,MMSEsimBer(:,1,1),'rx-','Linewidth',2);
grid on;hold on;
semilogy(Eb_No_dB,MMSEsimBer(:,1,2),'ro-','Linewidth',2);
semilogy(Eb_No_dB,MMSEsimBer(:,2,1),'bx-','Linewidth',2);
semilogy(Eb_No_dB,MMSEsimBer(:,2,2),'bo-','Linewidth',2);
semilogy(Eb_No_dB,MMSEsimBer(:,3,1),'kx-','Linewidth',2);
semilogy(Eb_No_dB,MMSEsimBer(:,3,2),'ko-','Linewidth',2);
set(gca,'FontSize',14);
hold off
axis square;
legend('9-tap MMSE For Channel 1', '15-tap MMSE For Channel 1','9-tap MMSE For Channel 2','15-tap MMSE For Channel 2',...
    '9-tap MMSE For Channel 3','15-tap MMSE For Channel 3','Location','southeast');
axis([0 40 0 0.6])
xlabel('Eb/No dB');
ylabel('Bit Error Rate');
title('BER For 3 Different Channels / MMSE Equalizer');