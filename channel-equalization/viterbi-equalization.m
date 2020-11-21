% Written in OCTAVE. Lines to comment or uncomment to make the code work in MATLAB were indicated.
% Below code implements Soft Decision Viterbi algorithm to equalize any given causal channel.
clear;clc;
bit_count = 1e3; % Number of bits in the stream for each SNR level
channel = [0.74 -0.514 0.37 0.216 0.062]; % Impulse response of the channel, it is assumed to be causal.
number_of_states = 2^(length(channel)-1); % Register has lenght(channel)-1 bits
binary_states = (dec2bin(0:number_of_states-1)); % Comment this line in MATLAB
binary_states = double(binary_states)-48;  % Comment this line in MATLAB
%binary_states = fliplr(de2bi(0:2^(length(channel)-1)-1)); % Uncomment this line in MATLAB
binary_states(binary_states==0) = -1;  % BPSK modulation
state_length = size(binary_states,2); % Register size
outputs = [[binary_states -ones(number_of_states,1)]*flipud(channel') [binary_states ones(number_of_states,1)]*flipud(channel')];
% For each state, two possible outputs for 0 or 1 inputs are stored in "outputs" matrix
SNR_start = 0;
SNR_increment = 2;
SNR_interval = SNR_start:SNR_increment:20;
err_store = zeros(length(SNR_interval),1); % Placeholder for errors
sim_limit = 10; % Number of simulations to average over
for iter = 1:length(SNR_interval)
  for sim_iter=1:sim_limit
      viterbi_metrics = zeros(number_of_states,1); % Placeholder to store metrics for each state
      viterbi_metrics(2:end) = Inf; % Stream starts from 1st state
      viterbi_path = zeros(number_of_states,bit_count); % Placeholder to store previous state
      noise_variance = 1/(2*10^(SNR_interval(iter)/10));
      data = 2.*randi([0 1],1,bit_count)-1;
      data =[-1.*ones(1,state_length) data];
      received = conv(data,channel); % Data with ISI
      received = received(1:length(data))+sqrt(noise_variance).*randn(1,length(data)); % Data with ISI and noise
      for bit_index=length(channel):bit_count+length(channel)-1
          state_metrics = abs(outputs-received(1,bit_index)); % Distances
          temp_viterbi_metrics = zeros(size(viterbi_metrics)); % Distance metrics for that bit_index
          for metric_index=1:number_of_states
              index1 = round(metric_index/2);  % 1st possible previous state 
              index2 = index1 + number_of_states/2; % 2nd possible previous state 
              index3 = mod(metric_index,2); % Possible input
              index3(index3 == 0) = 2; 
              alternative_values = [viterbi_metrics(index1)+state_metrics(index1,index3),viterbi_metrics(index2)+state_metrics(index2,index3)];
              [min_val, min_index] = min(alternative_values); % Pick the previous state that results in smaller total distance
              temp_viterbi_metrics(metric_index) = min_val;
              prev = [index1 index2];
              viterbi_path(metric_index,bit_index-length(channel)+1) = prev(min_index); % Store previous state
          end
          viterbi_metrics = temp_viterbi_metrics;
      end
      % Detection
      detected = zeros(1,bit_count);
      [exit_val, exit_index] = min(viterbi_metrics(:,end)); % Pick the state with the least total distance at the end.
      detected(end) = binary_states(exit_index,end);
      for detection_index = bit_count:-1:2
          temp = viterbi_path(exit_index,detection_index);
          detected(1,detection_index-1)=binary_states(temp,end);
          exit_index = temp;
      end
     err_store(iter) = err_store(iter)+length(find(detected~=data(length(channel):length(channel)+bit_count-1)));
   end 
   err_store(iter) = err_store(iter)/bit_count/sim_iter;
end
semilogy(SNR_interval,err_store,'*',"markersize",3);
title('Soft Decision Viterbi BER');
xlabel('E_b/N_o (dB)');
ylabel('BER');
grid on;
