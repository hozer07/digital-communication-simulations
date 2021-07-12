function [output_threshold1,output_threshold2,output_threshold3] = approx_threshold(tag_to_reader,source_to_tag,source_to_reader,tag_attenuation,N,noise_power)
mu = source_to_reader+ tag_attenuation.*source_to_tag.*tag_to_reader;
ss0_2 = 2/N.*abs(source_to_reader).^2*noise_power;
ss1_2 = 2/N.*abs(mu).^2*noise_power;
ssp_2 = ss0_2 + ss1_2;
gamma = abs(mu).^2 - abs(source_to_reader).^2;
output_threshold1 = abs(gamma)/2 + ssp_2/abs(gamma).*log(1+sqrt(1-exp(-gamma.^2./ssp_2)));
output_threshold2 = abs(gamma)/2;
p0 = @(x,ss0_2,ss1_2) 0.25./sqrt(pi.*ss0_2).*exp(-x.^2./4./ss0_2)+0.25./sqrt(pi.*ss1_2).*exp(-x.^2./4./ss1_2);
p1 = @(x,ss0_2,ss1_2,gamma) 1./sqrt(8.*pi.*(ss0_2+ss1_2)).*exp(-(x-gamma).^2./2./(ss0_2+ss1_2))+1./sqrt(8.*pi.*(ss0_2+ss1_2)).*exp(-(x+gamma).^2./2./(ss0_2+ss1_2));
lims = 0:1e-3:5;
p0_temp = p0(lims,ss0_2,ss1_2);
p1_temp = p1(lims,ss0_2,ss1_2,gamma);
temp = abs(p0_temp-p1_temp);
for j=2:length(lims)
    if(temp(j)>temp(j-1))
        break;
    end
end
output_threshold3 = lims(j);
end

