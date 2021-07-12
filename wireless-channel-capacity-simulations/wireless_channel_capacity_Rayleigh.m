clear;clc;
step_size = 1e-4;
SNR = 0:30;
lims = 0:step_size:8;
capacity_store = zeros(3,length(SNR));
for snr = 1:length(SNR)
    Es_No = 10^(SNR(snr)/10);
    gamma_i = Es_No.*lims;
    pdf = exp(-lims);
    capacity_store(1,snr) = sum(log2(1+gamma_i).*pdf).*step_size;
    capacity_store(2,snr) = log2(1+Es_No);
    gamma_zero = zeros(1,round(length(lims)/20));
    for i=1:round(length(lims)/20)
        gamma_zero(i) = sum((1/gamma_i(i)-1./gamma_i(i:end)).*pdf(i:end)).*step_size;
    end
    [~,index] = min(abs(1-gamma_zero));
    capacity_store(3,snr) = sum(log2(gamma_i(index:end)./gamma_i(index)).*pdf(index:end))*step_size;
end
plot(SNR,capacity_store(1,:),'-','LineWidth',2);grid on;hold on;xlabel('SNR (dB)');ylabel('Bits / Second / Hz');
plot(SNR,capacity_store(2,:),'-','LineWidth',1.5);
plot(SNR,capacity_store(3,:),'gd','MarkerSize',6);
legend('RX has Rayleigh CSI','AWGN Channel','TX & RX has Rayleigh CSI')
