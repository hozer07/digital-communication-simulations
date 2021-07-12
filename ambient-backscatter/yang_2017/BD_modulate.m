function [outputArg1] = BD_modulate(c,x,N,Nc,K,D_h)
outputArg1 = c;
for i=1:length(x)
    if(x(i)==1)
        indexes_first = D_h+(((i-1)*K*(N+Nc)+(N+Nc)/2+1):((i-1)*K*(N+Nc)+(N+Nc)));
        indexes = indexes_first;
        for iter = 2:K
            indexes = [indexes indexes_first+(iter-1)*(N+Nc)];
        end
%         if(indexes(end)>length(c))
%             crop = indexes(end)-length(c);
%             indexes(end-crop+1:end) = [];
%         end
        outputArg1(indexes) = -outputArg1(indexes);
    end
end