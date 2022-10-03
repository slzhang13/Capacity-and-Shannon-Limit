clear;
N = 16;
[x, w] = gaussHermite(N);
tic
for M = 8%[2 4 8 16 32 64]
    conste_norm = exp(1j*(0:M-1)*2*pi/M);
    K = log2(M);
    mapping = 0:M-1;%[0 1 3 2 6 7 5 4];
    labels = dec2bin(mapping, K)-'0';
    No = 1;
    SNR = -5:35;
    SNR_Lin = 10.^(0.1*SNR);
    C = zeros(size(SNR));
    
    for SNR_cnt = 1:length(SNR)

        snr_lin = SNR_Lin(SNR_cnt);
        Es = snr_lin;
        conste = sqrt(Es)*conste_norm;

        sum1 = 0;
        for k = 1:K
            for b = 0:1
                index = labels(:, k)==b;
                subset = conste(index);
                for cnt = 1:length(subset)
                    for nI = 1:N
                        for nQ = 1:N
                            fenzi = sum(exp(-abs(x(nI)+1j*x(nQ)+subset(cnt)-conste).^2));
                            fenmu = sum(exp(-abs(x(nI)+1j*x(nQ)+subset(cnt)-subset).^2));
                            sum1 = sum1+w(nI)*w(nQ)*log2(fenzi./fenmu);
                        end
                    end
                end
            end
            
        end
        C(SNR_cnt) = K-sum1/(pi*M);
    end
    plot(SNR,C,'LineWidth',2,'LineStyle','-.');
    grid on;
    hold on;
end
% thy_capacity = log2(1+SNR_Lin);
% plot(SNR,thy_capacity,'LineWidth',1);
toc
