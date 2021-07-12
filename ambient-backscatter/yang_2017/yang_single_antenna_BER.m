clear;clc;
N = 512;
Nc = 64;
c = 3*1e8;
fc = 2.4*1e9;
g_m = @(D_rd) c^2/(4*pi*fc^2*D_rd^2);
tau_g = 1; tau_f = 4; tau_h = 6; 
D_rd = 0.5; D_g = 0; D_f = 16; D_h = 16; D = 16; L = 22;
alpha = 0.3+0.4*1i;
symbols_per_frame = 20;
SNR_range = 0:5:30;
D_h_store = zeros(3,length(SNR_range));
D_store = zeros(3,length(SNR_range));
L_store = zeros(3,length(SNR_range));
threshold = @(J,SNR) (SNR+1)/(SNR.*(SNR+2)).*(SNR+sqrt(SNR.^2+(2.*SNR.*(SNR+2).*log(SNR+1))./J));
errs_K_1 = zeros(1,length(SNR_range));
theoretical_errs_K_1 = zeros(1,length(SNR_range)); 

K = 1;

for i=1:length(SNR_range)
    sim_count = 0;
    SNR = 10^(SNR_range(i)/10);
    noise_power = 2*abs(alpha)^2*g_m(D_rd)/SNR;
    legacy_data_syms = 2.*randi([0 1],1,K*N*symbols_per_frame)-1;
    legacy_data_with_prefix = zeros(1,K*(N+Nc)*symbols_per_frame);
    for j=1:K*symbols_per_frame
        legacy_data_syms((j-1)*N+1:j*N) = sqrt(N).*ifft(legacy_data_syms((j-1)*N+1:j*N));
        legacy_data_with_prefix((j-1)*(N+Nc)+1:(j-1)*(N+Nc)+Nc) = legacy_data_syms(j*N-Nc+1:j*N);
        legacy_data_with_prefix((j-1)*(N+Nc)+Nc+1:j*(N+Nc))=legacy_data_syms((j-1)*N+1:j*N);                 
    end
    J = Nc+D-L+1;
    while(errs_K_1(i)<1e3)      
        sim_count = sim_count + 1;
        source_to_BD_channel = channel_generator(tau_h);
        legacy_data_at_BD = conv([zeros(1,D_h) legacy_data_with_prefix],source_to_BD_channel);
        BD_data_bits = randi([0 1],1,symbols_per_frame);
        BD_sent = BD_modulate(legacy_data_at_BD,BD_data_bits,N,Nc,K,D_h);
        BD_to_receiver_channel = sqrt(0.5*g_m(D_rd)).*(randn+1i.*randn);
        detection_SNR = 2.*abs(alpha)^2.*abs(BD_to_receiver_channel).^2/noise_power;
        noise = sqrt(noise_power/2).*(randn(1,length(BD_sent))+1i.*randn(1,length(BD_sent)));
        BD_sent = alpha.*BD_to_receiver_channel.*BD_sent;
        source_to_receiver_channel = channel_generator(tau_f);
        source_to_receiver_signal = conv([zeros(1,D_f) legacy_data_with_prefix],source_to_receiver_channel);
        if(D_h+tau_h>D_f+tau_f)
            source_to_receiver_signal(end+1:end+D_h+tau_h-D_f-tau_f) = source_to_receiver_signal(end-(D_h+tau_h-D_f-tau_f)+1:end);
        elseif(D_h+tau_h<D_f+tau_f)
            BD_sent(end+1:end+D_f+tau_f-D_h-tau_h) = BD_sent(end-(D_f+tau_f-D_h-tau_h)+1:end);
        end
        total_received_signal = BD_sent + source_to_receiver_signal + noise;
        threshold_opt = threshold(K*J,detection_SNR);
        after_metric = metric_calculator(total_received_signal,N,Nc,symbols_per_frame,noise_power,D,L,K);
        detected = after_metric > threshold_opt;
        err = sum(detected~=BD_data_bits);
        errs_K_1(i)=errs_K_1(i)+err;
        theoretical_errs_K_1(i) = theoretical_errs_K_1(i)+(0.5*(qfunc(sqrt(K*J)*(threshold_opt-1))+qfunc(sqrt(K*J)*(1-threshold_opt/(detection_SNR+1)))));
    end
    errs_K_1(i) = errs_K_1(i)./sim_count./symbols_per_frame;
    theoretical_errs_K_1(i) = theoretical_errs_K_1(i)/sim_count;
end