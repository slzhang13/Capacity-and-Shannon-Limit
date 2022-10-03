clear;clc;

for M = [2 4 8 16]
    conste_norm = getConste(M);
    K = log2(M);
    No = 1;
    SNR = -10:0.5:30;
    SNR_Lin = 10.^(0.1*SNR);
    C = zeros(size(SNR));
    
    for SNR_cnt = 1:length(SNR)

        snr_lin = SNR_Lin(SNR_cnt);
        Es = snr_lin*0.5;
        conste = sqrt(Es)*conste_norm;
        sum1 = 0;
        t = sqrt(0.5)*randn(1,1e4);
        
        for i = 1:M
            
%             sum2 = 0;
%             for j = 1:M
%                 sum2 = sum2+exp(-(conste(j)-conste(i)).^2+...
%                             2*t.*(conste(j)-conste(i)));
%             end

            sum2 = sum(exp(-(conste.'-conste(i)).^2+...
                            2*t.*(conste.'-conste(i))),1);
            sum1 = sum1+mean(log2(sum2));
            
        end
        C(SNR_cnt) = K-sum1/M;
    end
    plot(SNR,C,'LineWidth',1,'LineStyle','-.');
    grid on;
    hold on;
end
thy_capacity = 0.5*log2(1+SNR_Lin);
plot(SNR,thy_capacity,'LineWidth',1);
