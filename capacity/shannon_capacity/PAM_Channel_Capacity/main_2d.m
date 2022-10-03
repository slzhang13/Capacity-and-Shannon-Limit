clear;clc;
N = 16;
[x, w] = gaussHermite(N);
tic
for M = [2 4 8 16 32 64]
    conste_norm = exp(1j*(0:M-1)*2*pi/M);
    K = log2(M);
    No = 1;
    SNR = -5:1:35;
    SNR_Lin = 10.^(0.1*SNR);
    C = zeros(size(SNR));
    
    for SNR_cnt = 1:length(SNR)

        snr_lin = SNR_Lin(SNR_cnt);
        Es = snr_lin;
        conste = sqrt(Es)*conste_norm;

        sum1 = 0;
        for i = 1:M
            sum2 = 0;
            for nI = 1:N
                for nQ = 1:N
                    sum3 = 0;
                    for j = 1:M
                        sum3 = sum3+exp(abs(x(nI)+1j*x(nQ))^2-abs(x(nI)+1j*x(nQ)+conste(i)-conste(j))^2);
                    end
                    sum2 = sum2+w(nI)*w(nQ)*log2(sum3);
                end
            end
            sum1 = sum1+sum2;
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
