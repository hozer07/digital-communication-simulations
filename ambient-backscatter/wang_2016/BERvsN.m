clear;clc;
K = 100;
tag_attenuation = 10^(-0.11);
SNR_range = (0:5:30);
num_sim = 1e5;
N_range = 10:10:200;
err_store_opt = zeros(1,length(N_range));
approx = zeros(1,length(N_range));
upper_bound = zeros(1,length(N_range));
lower_bound = zeros(1,length(N_range));
No = 10^(-3);
for n_var = 1:length(N_range)
    N = N_range(n_var);
    for sim = 1:num_sim
        source_to_tag = sqrt(0.5).*(randn+1i*randn);
        tag_to_reader = sqrt(0.5).*(randn+1i*randn);
        source_to_reader = sqrt(0.5).*(randn+1i*randn);
        mu_2 = abs(source_to_reader+tag_attenuation*tag_to_reader*source_to_tag)^2;
        h_2 = abs(source_to_reader)^2;
        low_var = min(mu_2,h_2);
        up_var = max(h_2,mu_2);
        delta_mu_h = abs(mu_2-h_2);
        var_mu_h = mu_2+h_2;
        upper_bound(n_var) = upper_bound(n_var)+ 1/8*exp(-delta_mu_h^2*N/No/32/up_var)+3/8*exp(-delta_mu_h^2*N/No/24/up_var);
        lower_bound(n_var) = lower_bound(n_var)+ 5/48*exp(-delta_mu_h^2*N/No/32/low_var)+5/16*exp(-delta_mu_h^2*N/No/24/low_var);
        approx(n_var) = approx(n_var) + 1/24*exp(-delta_mu_h^2*N/No/32/h_2)+1/8*exp(-delta_mu_h^2*N/No/24/h_2)+...
            1/24*exp(-delta_mu_h^2*N/No/32/mu_2)+1/8*exp(-delta_mu_h^2*N/No/24/mu_2)+1/24*exp(-delta_mu_h^2*N/No/16/var_mu_h)+...
            1/8*exp(-delta_mu_h^2*N/No/12/var_mu_h)-1/24*exp(-9*delta_mu_h^2*N/No/16/var_mu_h)-1/8*exp(-3*delta_mu_h^2*N/No/4/var_mu_h);
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
        [~,~,optimum_threshold] = approx_threshold(tag_to_reader,source_to_tag,source_to_reader,tag_attenuation,N,No);
        temp_detect = power_sum(y_n,K,N)';
        temp_detect = diff(temp_detect)';
        detected = (abs(temp_detect) >= optimum_threshold);
        err_store_opt(n_var) = err_store_opt(n_var)+length(find(detected~=a_k));
    end
end
err_store_opt = err_store_opt./num_sim./K;
upper_bound = upper_bound./num_sim;
lower_bound = lower_bound./num_sim;
approx = approx./num_sim;
semilogy(N_range,err_store_opt,'r-','LineWidth',1.5);hold on;
semilogy(N_range,upper_bound,'b-^','MarkerSize',5,'LineWidth',1.5);
semilogy(N_range,lower_bound,'b-v','MarkerSize',5,'LineWidth',1.5);
semilogy(N_range,approx,'g--o','MarkerSize',5,'LineWidth',1.5);
legend('Optimum threshold simulation','Upper bound','Lower bound','Approximate');
axis square;
grid on;ylim([1e-3 1e-1]);
set(gca,'FontSize',14);
xlabel('N');
ylabel('BER');