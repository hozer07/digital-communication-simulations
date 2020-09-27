clear;clc;
err_min=100;            % Minimum number of errors to proceed for a specified SNR level and modulation scheme
sim_lim = 20;           % Minimum number of simulations to average over for a specified SNR level and modulation scheme
bpsk_sym = [-1,1];      % Baseband BPSK symbols
qpsk_sym = pskmod(0:3,4,pi/4); % Baseband QPSK symbols
qpsk_bits =[0,0;0,1;1,1;1,0];  % Gray coded bits corresponding to each QPSK symbol
bfsk_sym = [1,1i];             % Baseband BFSK symbols 
SNR_int = -10:5:10;            % SNR levels
errs_store = zeros(3,length(SNR_int),2);     % Placeholder matrix to store errors for each SNR level and modulation scheme
ber = zeros(length(SNR_int),1);              % Placeholder array to store BER of QPSK
seq_len = 1e5;                 % Number of generated symbols for each simulation
for modulation =1:3            % i=1 --> QPSK, i=2 --> BPSK, i=3 --> BFSK
    for snr=SNR_int
        snr_index = snr/5+3;
        err_check = true;
        if(modulation==1)
            sim_num = 0;
            noise_variance = 1/(2*2*10.^(snr/10));      % Converted Eb/No to Es/No since Es = 1, then divided by 2 to have No/2
            while ((err_check==true || sim_num <= sim_lim))
                sim_num = sim_num + 1;
                info_seq = randi([0,3],seq_len,1);
                info_bits = qpsk_bits(info_seq+1,:);
                info_bits2 = info_bits';
                info_bits3 = info_bits2(:);
                four_sym_seq = conj(qpsk_sym(info_seq+1)');
                four_sym_seq = four_sym_seq + sqrt(noise_variance).*(randn(seq_len,1)+1i.*randn(seq_len,1));
                % Noise variance is divided by 2 to make overall noise
                % variance No/2, in-phase and quadrature components are
                % added
                four_sym_seq = repmat(four_sym_seq,[1,4]);
                norms = abs(four_sym_seq-repmat(qpsk_sym,[seq_len,1]));
                [~,detected] = min(norms,[],2); % Maximum Likelihood detection rule due to uniformly generated info bits
                detected = detected-1;
                err = length(find(detected~=info_seq));
                detected_bits = qpsk_bits(detected+1,:);
                detected_bits2 = detected_bits';
                detected_bits3 = detected_bits2(:);
                ber_temp = length(find(detected_bits3~=info_bits3));
                ber(snr_index) = ber(snr_index) + ber_temp;
                errs_store(modulation,snr_index,1)=errs_store(modulation,snr_index,1)+err;
                if(errs_store(modulation,snr_index,1)>100),err_check=false;end
            end
            errs_store(modulation,snr_index,2)=sim_num;
        else
            noise_variance = 1/(2*10.^(snr/10));
            err_check=true;
            sim_num=0;
            while (err_check==true || sim_num <= sim_lim)
                sim_num = sim_num + 1;
                info_seq = randi([0,1],seq_len,1);
                if(modulation==2)
                    sym_seq = bpsk_sym(info_seq+1)';
                    sym_seq = sym_seq + sqrt(noise_variance).*randn(seq_len,1);
                    detected = sym_seq>0;
                else
                    sym_seq = conj(bfsk_sym(info_seq+1)');
                    sym_seq = sym_seq+sqrt(noise_variance).*(randn(seq_len,1)+1i.*randn(seq_len,1));
                    sym_seq = repmat(sym_seq,[1,2]);
                    norms = abs(sym_seq-repmat(bfsk_sym,[seq_len,1]));
                    [~,detected] = min(norms,[],2);
                    detected = detected-1;                    
                end
                err = length(find(detected~=info_seq));
                errs_store(modulation,snr_index,1)=errs_store(modulation,snr_index,1)+err;
                if(errs_store(modulation,snr_index,1)>100),err_check=false;end           
            end
            errs_store(modulation,snr_index,2)=sim_num;
        end
    end
end
err_store = errs_store(:,:,1)./errs_store(:,:,2)./seq_len;
ber = ber'./errs_store(1,:,2)./seq_len./2;
EbNo = 10.^(SNR_int./10);
theoretical_symbol_error_qpsk = 2.*qfunc(sqrt(2.*EbNo)).*(1-0.5.*qfunc(sqrt(2.*EbNo)));
theoretical_bit_error_bpsk = qfunc(sqrt(2.*EbNo));
theoretical_bit_error_bfsk = qfunc(sqrt(EbNo));
semilogy(SNR_int,err_store(1,:),'r-s','MarkerSize',10);hold on;xlabel('Eb/No (dB)');ylabel('SER');grid on;
semilogy(SNR_int,err_store(2,:),'g-*','MarkerSize',10);title('Symbol Error Rate Plots');xticks(-10:5:10);
semilogy(SNR_int,err_store(3,:),'b');
semilogy(SNR_int,theoretical_symbol_error_qpsk,'r-x','MarkerSize',10);
semilogy(SNR_int,theoretical_bit_error_bpsk,'g-o','MarkerSize',10);
semilogy(SNR_int,theoretical_bit_error_bfsk,'b-p','MarkerSize',10);
legend('Simulation SER for QPSK','Simulation SER for BPSK','Simulation SER for BFSK',...
    'Theoretical SER for QPSK','Theoretical SER for BPSK','Theoretical SER for BFSK','Location','southwest');
set(gca,'FontSize',14);
figure,
semilogy(SNR_int,ber,'r-s','MarkerSize',10);hold on;xlabel('Eb/No (dB)');ylabel('BER');grid on;
semilogy(SNR_int,err_store(2,:),'g-*','MarkerSize',10);title('Bit Error Rate Plots');xticks(-10:5:10);
semilogy(SNR_int,err_store(3,:),'b-+');
semilogy(SNR_int,theoretical_symbol_error_qpsk./2,'r-x','MarkerSize',10);
semilogy(SNR_int,theoretical_bit_error_bpsk,'g-o','MarkerSize',10);
semilogy(SNR_int,theoretical_bit_error_bfsk,'b-p','MarkerSize',10);
legend('Simulation BER for QPSK','Simulation BER for BPSK','Simulation BER for BFSK',...
    'Theoretical BER for QPSK','Theoretical BER for BPSK','Theoretical BER for BFSK','Location','southwest');
set(gca,'FontSize',14);
