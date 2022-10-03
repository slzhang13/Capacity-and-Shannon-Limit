clear;
figure;
N = 16;
[x, w] = gaussHermite(N);
xx = x+1j*x.';
ww = w.*w.';

tic
for M = 256
%     conste_norm = exp(1j*(0:M-1)*2*pi/M);
conste_norm=qammod(0:M-1, M,'UnitAveragePower', true);
    K = log2(M);
    No = 1;
    SNR = -5:1:40;
    SNR_Lin = 10.^(0.1*SNR);
    C = zeros(size(SNR));
    
    for SNR_cnt = 1:length(SNR)

        snr_lin = SNR_Lin(SNR_cnt);
        Es = snr_lin;
        conste = sqrt(Es)*conste_norm;

        sum1 = 0;
        for i = 1:M
            sum2 = 0;
            for j = 1:M
                sum2 = sum2+exp(abs(xx).^2-abs(xx+conste(i)-conste(j)).^2);
            end
            sum2 = ww.*log2(sum2);
            sum1 = sum1+sum(sum(sum2));
        end
        C(SNR_cnt) = K-sum1/(pi*M);
    end
    plot(SNR,C,'LineWidth',1);
    grid on;
    hold on;
end

thy_capacity = log2(1+SNR_Lin);
plot(SNR,thy_capacity,'LineWidth',1);
toc