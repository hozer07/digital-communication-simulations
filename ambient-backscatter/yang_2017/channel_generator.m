function [channel] = channel_generator(delay_spread)
taps = randi([0 1],1,delay_spread-1);
channel = sqrt(0.5).*(randn(1,delay_spread)+1i.*randn(1,delay_spread)).*[1 taps].*exp(-2.*(0:delay_spread-1)./delay_spread);
channel = channel./sqrt(sum(abs(channel).^2));
end

