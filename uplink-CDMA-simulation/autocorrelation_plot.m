% Auto-Correlation Plot to demonstrate
% shift orthogonality of maximum length sequences
clear;clc;
[m_seq,~] = lfsr_generator(4,1);
m_seq(m_seq==1) = -1;
m_seq(m_seq==0) = 1;
corr_data = xcorr([zeros(1,2*length(m_seq)) m_seq zeros(1,2*length(m_seq))],repmat(m_seq,[1 5]));
plot(-74:74,corr_data./max(corr_data));title('Autocorrelation of m-sequence m=5');grid on;xlabel('nT_c');ylim([-0.5 1.2]);xlim([-35 35]);