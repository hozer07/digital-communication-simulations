clear;clc;
rayleigh_data = sqrt(0.5).*(randn(1,1e6)+1i.*randn(1,1e6));
rayleigh_power = abs(rayleigh_data).^2;
bin_width = 1e-2;
lims = 0:bin_width:4;
N = hist(abs(rayleigh_data),lims);hold on;
plot(lims,2.*lims.*exp(-lims.^2),'*','MarkerSize',8);
plot(lims,N./1e6./bin_width,'LineWidth',2);
legend('Rayleigh PDF','Envelope Histogram');
xlabel('z');ylabel('f_Z(z)');
axis square;
grid on;
bin_width = 1e-2;
lims2 = 0:bin_width:12;
N2 = hist(rayleigh_power,lims2);
figure,
hold on;
exp_pdf = exp(-lims2);
plot(lims2,exp_pdf,'*','MarkerSize',8);
plot(lims2(2:end),N2(2:end)./1e6./bin_width,'LineWidth',2.3);
legend('Exponential PDF','Power Histogram');
xlabel('v');ylabel('f_V(v)');
axis square;
grid on;