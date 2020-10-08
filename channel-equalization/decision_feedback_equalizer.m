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
DFEsimBer = zeros(length(Eb_No_dB),length(channels),length(tap_options)); % Placeholder for simulation results
detection = zeros(1,N); % Placeholder for detected symbols
for average = 1:5 % The number of simulations to average the results
    for snr=1:length(Eb_No_dB)
        infos = 2.*randi([0 1],1,N)-1; % Uniformly generated binary information sequence
        noise_variance = 0.5*(10^(-Eb_No_dB(snr)/10));
        noise = sqrt(noise_variance).*randn(1,N);
        for channel_index = 1:3
            channel_output = conv(infos,channels{channel_index},'same');
            y = channel_output+noise;
            system_length = length(channels{channel_index});
            for tap=1:length(tap_options)
                channel_matrix = toeplitz([channels{channel_index}(1);zeros(tap_options(tap)-1,1)],[channels{channel_index} zeros(1,tap_options(tap)-1)]);
                half = round(tap_options(tap)/2); % Feed forward filter lengths were chosen to be "half"
                ff_channel_matrix = toeplitz([channels{channel_index}(1);zeros(half-1,1)],[channels{channel_index}(1:ceil(system_length/2)) zeros(1,ceil(size(channel_matrix,2)/2)-ceil(system_length/2))]);
                w_ff = ((ff_channel_matrix*ff_channel_matrix'+ noise_variance*eye(half))\(ff_channel_matrix*[zeros(1,floor(size(channel_matrix,2)/2)) 1]'))';
                w_ff = [w_ff zeros(1,tap_options(tap)-length(w_ff))];
                y_out_df = conv(y,w_ff,'same');
                y_out_df = [y_out_df zeros(1,tap_options(tap)-half)];
                w_fb = channels{channel_index}(ceil(system_length/2)+1:end);
                w_fb = [w_fb zeros(1,half-1-length(w_fb))];
                w_fb = conv(w_fb,w_ff,'same');          
                for i=1:N
                    detected = (y_out_df(i)>0);
                    detection(1,i)=2*detected-1;
                    y_out_df(i+1:i+half-1) = y_out_df(i+1:i+half-1) - conv(w_fb,detection(1,i),'same');
                end
                DFEsimBer(snr,channel_index,tap)=DFEsimBer(snr,channel_index,tap)+sum((detection(1,:)~=infos));
            end
        end
    end
end
DFEsimBer = DFEsimBer./N./average; % simulated ber
figure
semilogy(Eb_No_dB,DFEsimBer(:,1,1),'rx-','Linewidth',2);
grid on;hold on;
semilogy(Eb_No_dB,DFEsimBer(:,1,2),'ro-','Linewidth',2);
semilogy(Eb_No_dB,DFEsimBer(:,2,1),'bx-','Linewidth',2);
semilogy(Eb_No_dB,DFEsimBer(:,2,2),'bo-','Linewidth',2);
semilogy(Eb_No_dB,DFEsimBer(:,3,1),'kx-','Linewidth',2);
semilogy(Eb_No_dB,DFEsimBer(:,3,2),'ko-','Linewidth',2);
set(gca,'FontSize',14);
hold off
axis square;
legend('9-tap DFE For Channel 1', '15-tap DFE For Channel 1','9-tap DFE For Channel 2','15-tap DFE For Channel 2',...
    '9-tap DFE For Channel 3','15-tap DFE For Channel 3','Location','southeast');
axis([0 40 0 0.6])
xlabel('Eb/No dB');
ylabel('Bit Error Rate');
title('BER For 3 Different Channels');