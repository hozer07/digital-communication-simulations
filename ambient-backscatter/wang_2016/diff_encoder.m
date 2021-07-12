function [outputArg1] = diff_encoder(inputArg1)
current = 0;
outputArg1 = zeros(1,length(inputArg1));
for i=1:length(inputArg1)
    outputArg1(i) = xor(current,inputArg1(i));
    current = outputArg1(i);
end
end

