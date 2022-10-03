SNR = -10:1:20;
SNR_Lin = 10.^(0.1*SNR);
tic
for M = 16
%     conste_tmp = [-3 -1 1 3]+1j*[3; 1; -1; -3];
% conste_tmp = conste_tmp(:);
% conste_norm = conste_tmp./sqrt(mean(abs(conste_tmp).^2));
conste_norm = qammod(0:M-1, M,'UnitAveragePower', true);
%     conste_tmp = [exp(1j*((0:4-1)*2*pi/4+pi/4)) 2.8*exp(1j*(2*pi*(0:12-1)/12+pi/12)) 5.3*exp(1j*(2*pi*(0:16-1)/16+pi/16))];
%     conste_norm = conste_tmp./sqrt(mean(abs(conste_tmp).^2));
    capacity = cal_capacity(conste_norm, 2, SNR);
    plot(SNR,capacity,'LineWidth',1,'LineStyle','-');
    grid on;
    hold on;
end
toc
thy_capacity = log2(1+SNR_Lin);
plot(SNR,thy_capacity,'LineWidth',1);