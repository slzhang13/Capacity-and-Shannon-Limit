clear;clc;
N = 16;
[x, w] = gaussHermite(N);
tic
for M = 8
    conste_norm = exp(1j*(0:M-1)*2*pi/M);
    K = log2(M);
    mapping = 0:M-1;
    labels = dec2bin(mapping, K)-'0';
    No = 1;
    SNR = -5:35;
    SNR_Lin = 10.^(0.1*SNR);
    C = zeros(size(SNR));
    
    for SNR_cnt = 1:length(SNR)

        snr_lin = SNR_Lin(SNR_cnt);
        Es = snr_lin;
        conste = sqrt(Es)*conste_norm;

        s = 0;
        for k = 1:K
            
            index_0 = labels(:, k)==0;
            index_1 = labels(:, k)==1;
            subset_0 = conste(index_0);
            subset_1 = conste(index_1);
            sum0 = 0;
            sum1 = 0;
            
            for cnt = 1:length(subset_0)
                
                for nI = 1:N
                    for nQ = 1:N
                        p0 = sum(exp(-abs(x(nI)+1j*x(nQ)+subset_0(cnt)-subset_0).^2));
                        p1 = sum(exp(-abs(x(nI)+1j*x(nQ)+subset_0(cnt)-subset_1).^2));
                        p = p0+p1;
                        sum0 = sum0+w(nI)*w(nQ)*log2(p/p0);
                        
                        p0_ = sum(exp(-abs(x(nI)+1j*x(nQ)+subset_1(cnt)-subset_0).^2));
                        p1_ = sum(exp(-abs(x(nI)+1j*x(nQ)+subset_1(cnt)-subset_1).^2));
                        p_ = p0_+p1_;
                        sum1 = sum1+w(nI)*w(nQ)*log2(p_/p1_);
                    end
                end
                
            end
            s = s+(sum0+sum1)/(M*pi);
            
        end
        C(SNR_cnt) = K-s;
    end
    
    plot(SNR,C,'LineWidth',2,'LineStyle','-.');
    grid on;
    hold on;
end
% thy_capacity = log2(1+SNR_Lin);
% plot(SNR,thy_capacity,'LineWidth',1);
toc
