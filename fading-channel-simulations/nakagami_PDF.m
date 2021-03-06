clear;clc;
% K=0.5;
rician_data = sqrt(0.5)*0.6+1i*sqrt(0.5)*0.8+sqrt(0.5).*(randn(1,1e6)+1i.*(randn(1,1e6)));
rician_power = abs(rician_data).^2;
bin_width = 1e-2;
lims = 0:bin_width:4;
m = @(K) (K+1).^2/(2*K+1);
K=0.5;
N = hist(abs(rician_data),lims);hold on;
plot(lims,2*m(K).^m(K).*lims.^(2*m(K)-1).*exp(-m(K).*lims.^2./(K+1))./(gamma(m(K)).*(K+1)^(m(K))),'*','MarkerSize',8);
plot(lims,N./1e6./bin_width,'LineWidth',2);
legend('Nakagami Envelope PDF K=0.5','Envelope Histogram');
xlabel('z');ylabel('f_Z(z)');
axis square;
grid on;
figure,
lims = 0:bin_width:12;
N = hist(rician_power,lims);hold on;
plot(lims,(m(K)/(K+1)).^m(K).*lims.^(m(K)-1)./gamma(m(K)).*exp(-m(K).*lims./(K+1)),'*','MarkerSize',8);
plot(lims,N./1e6./bin_width,'LineWidth',2);
legend('Nakagami Power PDF K=0.5','Power Histogram');
xlabel('z');ylabel('f_V(v)');
axis square;
grid on;
% K = 1;
lims = 0:bin_width:6;
figure,
rician_data = 0.6+1i*0.8+sqrt(0.5).*(randn(1,1e6)+1i.*(randn(1,1e6)));
rician_power = abs(rician_data).^2;
K = 1;
N = hist(abs(rician_data),lims);
hold on;
plot(lims,2*m(K).^m(K).*lims.^(2*m(K)-1).*exp(-m(K).*lims.^2./(K+1))./(gamma(m(K)).*(K+1)^(m(K))),'*','MarkerSize',8);
plot(lims,N./1e6./bin_width,'LineWidth',2);
legend('Nakagami Envelope PDF K=1','Envelope Histogram');
xlabel('z');ylabel('f_Z(z)');
axis square;
grid on;
figure,
lims = 0:bin_width:15;
N = hist(rician_power,lims);hold on;
plot(lims,(m(K)/(K+1)).^m(K).*lims.^(m(K)-1)./gamma(m(K)).*exp(-m(K).*lims./(K+1)),'*','MarkerSize',8);
plot(lims,N./1e6./bin_width,'LineWidth',2);
legend('Nakagami Power PDF K=1','Power Histogram');
xlabel('z');ylabel('f_V(v)');
axis square;
grid on;
figure,
% K = 3;
rician_data = 0.6*sqrt(3)+1i*0.8*sqrt(3)+sqrt(0.5).*(randn(1,1e6)+1i.*(randn(1,1e6)));
rician_power = abs(rician_data).^2;
lims = 0:bin_width:6;
N = hist(abs(rician_data),lims);
K = 3;
hold on;
plot(lims,2*m(K).^m(K).*lims.^(2*m(K)-1).*exp(-m(K).*lims.^2./(K+1))./(gamma(m(K)).*(K+1)^(m(K))),'*','MarkerSize',8);
plot(lims,N./1e6./bin_width,'LineWidth',2);
legend('Nakagami Envelope PDF K=3','Envelope Histogram');
xlabel('z');ylabel('f_Z(z)');
axis square;
grid on;
figure,
lims = 0:bin_width:25;
N = hist(rician_power,lims);hold on;
plot(lims,(m(K)/(K+1)).^m(K).*lims.^(m(K)-1)./gamma(m(K)).*exp(-m(K).*lims./(K+1)),'*','MarkerSize',8);
plot(lims,N./1e6./bin_width,'LineWidth',2);
legend('Nakagami Power PDF K=3','Power Histogram');
xlabel('z');ylabel('f_V(v)');
axis square;
grid on;
figure,
% K =10;
lims = 0:bin_width:10;
rician_data = 0.6*sqrt(10)+1i*sqrt(10)*0.8+sqrt(0.5).*(randn(1,1e6)+1i.*(randn(1,1e6)));
rician_power = abs(rician_data).^2;
N = hist(abs(rician_data),lims);
hold on;
K=10;
plot(lims,2*m(K).^m(K).*lims.^(2*m(K)-1).*exp(-m(K).*lims.^2./(K+1))./(gamma(m(K)).*(K+1)^(m(K))),'*','MarkerSize',8);
plot(lims,N./1e6./bin_width,'LineWidth',2);
legend('Nakagami Envelope PDF K=10','Envelope Histogram');
xlabel('z');ylabel('f_Z(z)');
axis square;
grid on;
figure,
lims = 0:bin_width:50;
N = hist(rician_power,lims);hold on;
plot(lims,(m(K)/(K+1)).^m(K).*lims.^(m(K)-1)./gamma(m(K)).*exp(-m(K).*lims./(K+1)),'*','MarkerSize',8);
plot(lims,N./1e6./bin_width,'LineWidth',2);
legend('Nakagami Power PDF K=10','Power Histogram');
xlabel('z');ylabel('f_V(v)');
axis square;
grid on;