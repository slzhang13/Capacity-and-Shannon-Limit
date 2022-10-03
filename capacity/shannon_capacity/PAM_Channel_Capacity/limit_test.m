clear;
N = 16;
[x, w] = gaussHermite(N);
xx = x+1j*x.';
ww = w.*w.';
M = 16;
conste_tmp = [-3 -1 1 3]+1j*[3; 1; -1; -3];
conste_tmp = conste_tmp(:);
% conste_norm = exp(1j*(0:M-1)*2*pi/M);
% conste_tmp = [0.186*exp(1j*((0:4-1)*2*pi/4+pi/4)) 0.382*exp(1j*(2*pi*(0:8-1)/8+pi/8)) 0.560*exp(1j*(2*pi*(0:8-1)/8)) 0.754*exp(1j*(2*pi*(0:12-1)/12+pi/12)) 1.034*exp(1j*(2*pi*(0:16-1)/16+pi/16)) 1.506*exp(1j*(2*pi*(0:16-1)/16))];
conste_norm = conste_tmp./sqrt(mean(abs(conste_tmp).^2));
% conste_norm = qammod(0:M-1, M,'UnitAveragePower', true);
K = log2(M);
No = 1;
R_vec = 2/3;
limits = zeros(size(R_vec));
for R_cnt = 1:length(R_vec)
    R = R_vec(R_cnt);
    snr_start = -20;
    snr_end = 40;
    snr_num = 100;
    snr_dB = linspace(snr_start,snr_end,snr_num);
%     snr_Lin = 10.^(0.1*snr_dB);
%     C = zeros(size(snr_dB));
    iter = 1;
    while(1)
        fprintf("%d\n",iter);
        iter = iter+1;
        C = cal_capacity(conste_norm, 2, snr_dB);
        %[3 11 9 1 7 6 2 10 14 15 13 12 8 0 4 5]
%         for cnt = 1:length(snr_dB)
% 
%             snr = snr_Lin(cnt);
%             Es = snr;
%             conste = sqrt(Es)*conste_norm;
% 
%             sum1 = 0;
%             for i = 1:M
%                 sum2 = 0;
%                 for j = 1:M
%                     sum2 = sum2+exp(abs(xx).^2-abs(xx+conste(i)-conste(j)).^2);
%                 end
%                 sum2 = ww.*log2(sum2);
%                 sum1 = sum1+sum(sum(sum2));
%             end
%             C(cnt) = K-sum1/(pi*M);
%         end

        if length(find(C==R*K))==1
            fprintf("limit=%d\n",snr_dB(C==R*K));
            limits(R_cnt) = snr_dB(C==R*K);
            break;
        end

        tmp = find(C>R*K);
        g_index = tmp(1);
        l_index = tmp(1)-1;
        if(C(g_index)-C(l_index)<=1e-6)
            fprintf("limit=%d\n",snr_dB(l_index));
            limits(R_cnt) = snr_dB(l_index);
            break;
        end
        snr_start = snr_dB(l_index);
        snr_end = snr_dB(g_index);
        snr_dB = linspace(snr_start,snr_end,snr_num);
        snr_Lin = 10.^(0.1*snr_dB);
        C = zeros(size(snr_dB));
    end
end