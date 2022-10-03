clear;
N = 16;
[x, w] = gaussHermite(N);
tic
p = 2.^[-6 -6 -5 -4 -4 -4 -3 -3 -3 -3 -4 -4 -4 -5 -6 -6];
% p = ones(1,16)/16;
for M = 16
    conste_norm = getConste(M);
    mean_Es = sum(conste_norm.^2.*p);
    conste_norm = conste_norm./sqrt(mean_Es);
    K = log2(M);
    No = 1;
    SNR = -10:1:100;
    SNR_Lin = 10.^(0.1*SNR);
    C = zeros(size(SNR));
    
    for SNR_cnt = 1:length(SNR)

        snr_lin = SNR_Lin(SNR_cnt);
        Es = snr_lin*0.5;
        conste = sqrt(Es)*conste_norm;

        sum1 = 0;
        for i = 1:M
            sum2 = 0;
            for n = 1:N
                sum3 = 0;
                for j = 1:M
                    sum3 = sum3+p(j)*exp(-(conste(j)-conste(i)).^2+...
                            2*x(n)*(conste(j)-conste(i)));
                end
                sum2 = sum2+w(n)*log2(sum3);
            end
            sum1 = sum1+sum2*p(i);
        end
        C(SNR_cnt) = -sum1/(sqrt(pi));
    end
    plot(SNR,C,'LineWidth',1);
    grid on;
    hold on;
end
% thy_capacity = 0.5*log2(1+SNR_Lin);
% plot(SNR,thy_capacity,'LineWidth',1);
toc