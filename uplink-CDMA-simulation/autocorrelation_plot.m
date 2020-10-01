% Auto-Correlation Plot to demonstrate
% shift orthogonality of maximum length sequences
clear;clc;
[m_seq,~] = lfsr_generator(5,1);
m_seq(m_seq==1) = -1;
m_seq(m_seq==0) = 1;
corr_data = xcorr([zeros(1,3*length(m_seq)) m_seq zeros(1,3*length(m_seq))],repmat(m_seq,[1 7]));
plot(-floor(length(corr_data)/2):floor(length(corr_data)/2),corr_data./max(corr_data));
title('Autocorrelation of m-sequence m=4');grid on;xlabel('nT_c');
ylim([-0.5 1.2]);
xlim([-35 35]);
