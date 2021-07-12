clear;clc;
rician_data = sqrt(0.5)*0.6+1i*sqrt(0.5)*0.8+sqrt(0.5).*(randn(1,1e6)+1i.*(randn(1,1e6)));
bin_width = 1e-2;
lims = 0:bin_width:4;
N = hist(abs(rician_data),lims);hold on;
plot(lims,2.*lims.*exp(-(lims.^2+0.5)).*besseli(0,lims.*2.*sqrt(0.5)),'*','MarkerSize',8);
plot(lims,N./1e6./bin_width,'LineWidth',2);
plot(lims,2.*lims.*exp(-lims.^2),'.','MarkerSize',7);
legend('Rician PDF K=0.5','Envelope Histogram','Rayleigh PDF');
xlabel('z');ylabel('f_Z(z)');
axis square;
grid on;
figure,
rician_data = 0.6+1i*0.8+sqrt(0.5).*(randn(1,1e6)+1i.*(randn(1,1e6)));
N = hist(abs(rician_data),lims);
hold on;
plot(lims,2.*lims.*exp(-(lims.^2+1)).*besseli(0,lims.*2),'*','MarkerSize',8);
plot(lims,N./1e6./bin_width,'LineWidth',2);
plot(lims,2.*lims.*exp(-lims.^2),'.','MarkerSize',7);
legend('Rician PDF K=1','Envelope Histogram','Rayleigh PDF');
xlabel('z');ylabel('f_Z(z)');
axis square;
grid on;
figure,
rician_data = 0.6*sqrt(3)+1i*0.8*sqrt(3)+sqrt(0.5).*(randn(1,1e6)+1i.*(randn(1,1e6)));
lims = 0:bin_width:6;
N = hist(abs(rician_data),lims);
hold on;
plot(lims,2.*lims.*exp(-(lims.^2+3)).*besseli(0,lims.*2.*sqrt(3)),'*','MarkerSize',8);
plot(lims,N./1e6./bin_width,'LineWidth',2);
plot(lims,2.*lims.*exp(-lims.^2),'.','MarkerSize',7);
legend('Rician PDF K=3','Envelope Histogram','Rayleigh PDF');
xlabel('z');ylabel('f_Z(z)');
axis square;
grid on;
figure,
lims = 0:bin_width:7;
rician_data = 0.6*sqrt(10)+1i*sqrt(10)*0.8+sqrt(0.5).*(randn(1,1e6)+1i.*(randn(1,1e6)));
N = hist(abs(rician_data),lims);
hold on;
plot(lims,2.*lims.*exp(-(lims.^2+10)).*besseli(0,lims.*2.*sqrt(10)),'*','MarkerSize',8);
plot(lims,N./1e6./bin_width,'LineWidth',2);
plot(lims,2.*lims.*exp(-lims.^2),'.','MarkerSize',7);
legend('Rician PDF K=10','Envelope Histogram','Rayleigh PDF');
xlabel('z');ylabel('f_Z(z)');
axis square;
grid on;