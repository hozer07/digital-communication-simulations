clear;clc;
K = 100;
tag_attenuation = 10^(-0.11);
SNR_range = (0:5:30);
num_sim = 1e4;
err_store_exp = zeros(2,length(SNR_range));
err_store_gamma_over_2 = zeros(2,length(SNR_range));
err_store_appr = zeros(2,length(SNR_range));
err_store_opt = zeros(2,length(SNR_range));
thresholds_opt = zeros(1,length(SNR_range));
thresholds_exp = zeros(1,length(SNR_range));
thresholds_appr = zeros(1,length(SNR_range));
thresholds_gamm = zeros(1,length(SNR_range));
thresholds_eq = zeros(1,length(SNR_range));
N_range = [20 50];
for snr=1:length(SNR_range)
    No = 10^(-SNR_range(snr)/10)/2;
    for n_var = [1,2]
        N = N_range(n_var);
        for sim = 1:num_sim
            source_to_tag = sqrt(0.5).*(randn+1i*randn);
            tag_to_reader = sqrt(0.5).*(randn+1i*randn);
            source_to_reader = sqrt(0.5).*(randn+1i*randn);
            s_n = exp(1i*pi/2.*randi([0 3],1,K*N)+1i*pi/4);
            x_n = source_to_tag.*s_n;
            a_k = randi([0 1],1,K);
            b_k = diff_encoder(a_k);
            b_n = repelem(b_k,1,N);
            a_n = tag_attenuation.*b_n.*x_n;
            noise = sqrt(No/2).*(randn(1,K*N)+1i.*randn(1,K*N));
            y_n = source_to_reader.*s_n+tag_to_reader.*a_n+noise;
            initial_dummy = sqrt(No/2).*(randn(1,N)+1i.*randn(1,N)) + source_to_reader.*exp(1i*pi/2.*randi([0 3],1,N)+1i*pi/4);
            y_n = [initial_dummy y_n];
            [approximate_threshold,gamma_over_2,optimum_threshold] = approx_threshold(tag_to_reader,source_to_tag,source_to_reader,tag_attenuation,N,No);
            temp_detect = power_sum(y_n,K,N)';
            temp_detect = diff(temp_detect)';
            exp_threshold = mean(abs(temp_detect)); 
            detected = (abs(temp_detect) >= exp_threshold);
            err_store_exp(n_var,snr) = err_store_exp(n_var,snr)+length(find(detected~=a_k));
            detected = (abs(temp_detect) >= approximate_threshold);
            err_store_appr(n_var,snr) = err_store_appr(n_var,snr)+length(find(detected~=a_k));
            detected = (abs(temp_detect) >= gamma_over_2);
            err_store_gamma_over_2(n_var,snr) = err_store_gamma_over_2(n_var,snr)+length(find(detected~=a_k));
            detected = (abs(temp_detect) >= optimum_threshold);
            err_store_opt(n_var,snr) = err_store_opt(n_var,snr)+length(find(detected~=a_k));
        end
    end
end
err_store_exp = err_store_exp./num_sim./K;
err_store_appr = err_store_appr./num_sim./K;
err_store_gamma_over_2 = err_store_gamma_over_2./num_sim./K;
err_store_opt = err_store_opt./num_sim./K;
semilogy(SNR_range,err_store_exp(1,:),'b-s','MarkerSize',10,'LineWidth',2);hold on;grid on;
semilogy(SNR_range,err_store_exp(2,:),'b--*','MarkerSize',10,'LineWidth',2);
semilogy(SNR_range,err_store_appr(1,:),'r-d','MarkerSize',7,'LineWidth',2);
semilogy(SNR_range,err_store_appr(2,:),'r--v','MarkerSize',7,'LineWidth',2);
semilogy(SNR_range,err_store_gamma_over_2(1,:),'k-+','MarkerSize',5,'LineWidth',1.5);
semilogy(SNR_range,err_store_gamma_over_2(2,:),'k--x','MarkerSize',5,'LineWidth',1.5);
semilogy(SNR_range,err_store_opt(1,:),'g-o','MarkerSize',5,'LineWidth',1.5);
semilogy(SNR_range,err_store_opt(2,:),'g--o','MarkerSize',5,'LineWidth',1.5);
legend('N=20, Th=E[|\Phi|]','N=50, Th=E[|\Phi|]','N=20, Th=approx','N=50, Th=approx','N=20, Th=|\delta|/2','N=50, Th=|\delta|/2','N=20, Th=opt','N=50, Th=opt');
axis square;
set(gca,'FontSize',14);
xlabel('SNR (dB)');
ylabel('BER');
ylim([1e-3 1]);