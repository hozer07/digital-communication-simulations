function [R] = metric_calculator(y_n,N,Nc,symbols_per_frame,noise_power,D,L,K)
z = y_n - [y_n(1+N:end) zeros(1,N)];
R = zeros(1,symbols_per_frame);
J = Nc+D-L+1;
for i=1:symbols_per_frame
    start_index = (L:Nc+D)+(i-1)*K*(N+Nc);
    indexes = start_index;
    for K_index = 2:K
        indexes = [indexes start_index+(K_index-1)*(N+Nc)];
    end
    R(i) = sum(abs(z(indexes)).^2)./(K*J*noise_power*2);
end
end

