clear;clc;
K = 100;
tag_attenuation = 10^(-0.11);
SNR_range = (0:5:30);
num_sim = 1e4;
err_store_opt_0_1 = zeros(1,length(SNR_range));
err_store_opt_1_0 = zeros(1,length(SNR_range));
err_store_eq_0_1 = zeros(1,length(SNR_range));
err_store_eq_1_0 = zeros(1,length(SNR_range));
Be = @ (ss0_2,ss1_2,gamma) ((2+sqrt(2)).*0.717.*sqrt(max(ss0_2,ss1_2))+4*0.416*abs(gamma));
Ce = @ (ss0_2,ss1_2,gamma) (2*0.416*gamma^2+2*max(ss0_2,ss1_2)*log(2)+2*0.717*abs(gamma)*sqrt(max(ss0_2,ss1_2)));
for snr=1:length(SNR_range)
    No = 10^(-SNR_range(snr)/10);
    N = 20;
    for sim = 1:num_sim
        source_to_tag = sqrt(0.5).*(randn+1i*randn);
        tag_to_reader = sqrt(0.5).*(randn+1i*randn);
        source_to_reader = sqrt(0.5).*(randn+1i*randn);
        [~,~,optimum_threshold] = approx_threshold(tag_to_reader,source_to_tag,source_to_reader,tag_attenuation,N,No);
        mu = source_to_reader+tag_attenuation*tag_to_reader*source_to_tag;
        h = source_to_reader;
        ss0_2 = 2/N.*abs(source_to_reader).^2*No;
        ss1_2 = 2/N.*abs(mu).^2*No;
        gamma = abs(mu).^2 - abs(source_to_reader).^2;
        if(abs(ss0_2-ss1_2)>=1)
            Th_eq = 1/0.832*(Be(ss0_2,ss1_2,gamma)-sqrt(Be(ss0_2,ss1_2,gamma)^2-4*0.416*Ce(ss0_2,ss1_2,gamma)));
        else
            Th_eq = abs(gamma)/2 + log(2)*(ss0_2+ss1_2)/(2*0.717*sqrt(ss0_2+ss1_2)+2*0.416*abs(gamma));
        end
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
        temp_detect = power_sum(y_n,K,N)';
        temp_detect = diff(temp_detect)';
        detected = (abs(temp_detect) >= optimum_threshold);
        err_temp = (detected~=a_k & detected==1);
        err_store_opt_0_1(snr) = err_store_opt_0_1(snr) + sum(err_temp);
        err_temp = (detected~=a_k & detected==0);
        err_store_opt_1_0(snr) = err_store_opt_1_0(snr) + sum(err_temp);
        detected = (abs(temp_detect) >= Th_eq);
        err_temp = (detected~=a_k & detected==1);
        err_store_eq_0_1(snr) = err_store_eq_0_1(snr) + sum(err_temp);
        err_temp = (detected~=a_k & detected==0);
        err_store_eq_1_0(snr) = err_store_eq_1_0(snr) + sum(err_temp);
    end
end
err_store_opt_0_1 = err_store_opt_0_1./num_sim./K;
err_store_opt_1_0 = err_store_opt_1_0./num_sim./K;
err_store_eq_0_1 = err_store_eq_0_1./num_sim./K;
err_store_eq_1_0 = err_store_eq_1_0./num_sim./K;
semilogy(SNR_range,err_store_opt_0_1,'r--x','MarkerSize',7,'LineWidth',2);hold on;grid on;
semilogy(SNR_range,err_store_opt_1_0,'r--+','MarkerSize',7,'LineWidth',2);
semilogy(SNR_range,err_store_eq_0_1,'b--d','MarkerSize',7,'LineWidth',2);
semilogy(SNR_range,err_store_eq_1_0,'b--v','MarkerSize',7,'LineWidth',2);
semilogy(SNR_range,err_store_eq_0_1+err_store_eq_1_0,'b-s','MarkerSize',9,'LineWidth',2);
semilogy(SNR_range,err_store_opt_0_1+err_store_opt_1_0,'r-o','MarkerSize',5,'LineWidth',2);
legend('opt 0->1','opt 1->0','Th_e_q 0->1','Th_e_q 1->0','eq 0+1','opt 0+1');
axis square;
set(gca,'FontSize',14);
xlabel('SNR (dB)');
ylabel('BER');
ylim([1e-3 1]);