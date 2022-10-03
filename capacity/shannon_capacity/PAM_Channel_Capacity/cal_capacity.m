function capacity = cal_capacity(conste_norm, dim, SNR, labels)

    M = length(conste_norm);
    K = log2(M);
    No = 1;
    SNR_Lin = 10.^(0.1*SNR);
    capacity = zeros(size(SNR));

    N = 16;
    [x, w] = gaussHermite(N);
    xx = x+1j*x.';
    ww = w.*w.';

    if (nargin==4)

        bin_labels = dec2bin(labels, K)-'0';
        for SNR_cnt = 1:length(SNR)

            snr_lin = SNR_Lin(SNR_cnt);
            Es = snr_lin*No*dim*0.5;
            conste = sqrt(Es)*conste_norm;
            s = 0;
            for k = 1:K
                
                index_0 = bin_labels(:, k)==0;
                index_1 = bin_labels(:, k)==1;
                subset_0 = reshape(conste(index_0),1,1,M/2);
                subset_1 = reshape(conste(index_1),1,1,M/2);
                sum0 = 0;
                sum1 = 0;
                
                for cnt = 1:length(subset_0)
                    
                    p0 = sum(exp(-abs(xx+subset_0(cnt)-subset_0).^2), 3);
                    p1 = sum(exp(-abs(xx+subset_0(cnt)-subset_1).^2), 3);
                    p = p0+p1;
                    sum0 = sum0+ww.*log2(p./p0);

                    p0_ = sum(exp(-abs(xx+subset_1(cnt)-subset_0).^2), 3);
                    p1_ = sum(exp(-abs(xx+subset_1(cnt)-subset_1).^2), 3);
                    p_ = p0_+p1_;
                    sum1 = sum1+ww.*log2(p_./p1_);
                    
                end
                tmp = sum0+sum1;
                s = s+sum(tmp(:))/(M*pi);
                
            end
            capacity(SNR_cnt) = K-s;
        end
    else
        for SNR_cnt = 1:length(SNR)

            snr_lin = SNR_Lin(SNR_cnt);
            Es = snr_lin*No*dim*0.5;
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
            capacity(SNR_cnt) = K-sum1/(pi*M);
        end
    end

function [x, w] = gaussHermite(n)

    p = zeros(n+1);
    p(1, end) = 1;
    p(2, n:end) = [1 0]*2;
    for k = 3:n+1
        p(k, n-k+2:end) = 2*[p(k-1, n-k+3:end) 0]-2*(k-2)*[0 0 p(k-2, n-k+4:end)];
    end

    x = roots(p(n+1, :));
    w = zeros(1, n);
    for i = 1:n
        w(i) = 2^(n-1)*factorial(n)*sqrt(pi)/(n^2*polyval(p(n, :), x(i))^2);
    end
    x = x.';
    