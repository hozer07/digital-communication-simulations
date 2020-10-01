function [outputArg,outputArg2] = lfsr_generator(m,users)
    if(m==4)
        poly = [1 0 0 1];
        seed = [1 0 0 0];
        A = [poly;[eye(3) [0;0;0]]];
    elseif(m==5)
        poly = [0 1 0 0 1];
        A = [poly;[eye(4) [0;0;0;0]]];
        seed = [0 0 0 0 1];      
    end
    sequences = zeros(2^m,m);
    sequences(1,:) = seed;
    seq = sequences(1,end);
    for i=2:2^m
        sequences(i,:) = mod(A*sequences(i-1,:)',2);
        seq = [seq sequences(i,end)];
    end
    outputArg = seq(1:end-1);
    outputArg2 = outputArg;
    if(users>1)
        for i=2:users
            outputArg2(i,:) = [outputArg2(i-1,end) outputArg2(i-1,1:end-1)];
        end
    end
end