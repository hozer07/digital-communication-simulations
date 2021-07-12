function [outputArg1] = power_sum(inputArg1,K,N)
    temp = reshape(inputArg1,N,K+1);
    outputArg1 = sum(abs(temp).^2,1)./N;   
end

