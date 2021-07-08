clear;clc;
step_size = 1e-4;
SNR = 0:30;
lims = 0:step_size:8;
capacity_store = zeros(3,length(SNR));
m = 2;
for snr = 1:length(SNR)
    Es_No = 10^(SNR(snr)/10);
    pdf = (m).^m.*lims.^(m-1)./gamma(m).*exp(-m.*lims);
    gamma_i = Es_No.*lims;
    capacity_store(1,snr) = sum(log2(1+gamma_i).*pdf).*step_size;
    capacity_store(2,snr) = log2(1+Es_No);
    integral_result = zeros(1,round(length(lims)/20));
    for i=1:round(length(lims)/20)
        integral_result(i) = sum((1/gamma_i(i)-1./gamma_i(i:end)).*pdf(i:end)).*step_size;
    end
    [~,index] = min(abs(1-integral_result));
    capacity_store(3,snr) = sum(log2(gamma_i(index:end)./gamma_i(index)).*pdf(index:end))*step_size;
end
plot(SNR,capacity_store(1,:),'-','LineWidth',2);grid on;hold on;xlabel('SNR (dB)');ylabel('Bits / Second / Hz');
plot(SNR,capacity_store(2,:),'-','LineWidth',1.5);
plot(SNR,capacity_store(3,:),'gd','MarkerSize',6);
legend('RX has Nakagami m=2 CSI','AWGN Channel','TX & RX has Nakagami m=2 CSI')