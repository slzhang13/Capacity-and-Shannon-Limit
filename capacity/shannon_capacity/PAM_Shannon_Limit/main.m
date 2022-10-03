clear;clc;
N = 16;
[x, w] = gaussHermite(N);
tic
M = 16;
conste_norm = getConste(M);
K = log2(M);
No = 1;
R_vec = 1/2;%0.01:0.01:0.99;
limits = zeros(size(R_vec));
for R_cnt = 1:length(R_vec)
R = R_vec(R_cnt);
ebno_start = -2;
ebno_end = 20;
ebno_num = 100;
EbNo_dB = linspace(ebno_start,ebno_end,ebno_num);
EbNo_Lin = 10.^(0.1*EbNo_dB);
C = zeros(size(EbNo_dB));
iter = 1;
p = 2.^[-6 -6 -5 -4 -4 -4 -3 -3 -3 -3 -4 -4 -4 -5 -6 -6];
while(1)
    fprintf("%d\n",iter);
    iter = iter+1;
    
    for cnt = 1:length(EbNo_dB)

        Eb = EbNo_Lin(cnt)*No;
        Es = Eb*K*R;
        conste = sqrt(Es)*conste_norm;

        sum1 = 0;
        for i = 1:M
            sum2 = 0;
            for n = 1:N
                sum3 = 0;
                for j = 1:M
                    sum3 = sum3+exp(-(conste(j)-conste(i)).^2+...
                            2*x(n)*(conste(j)-conste(i)));
                end
%                 sum2 = sum2+w(n)*log2(sum3);
                sum2 = sum2+w(n)*log2(sum3)*p(i);
            end
            sum1 = sum1+sum2;
        end
%         C(cnt) = K-sum1/(sqrt(pi)*M);
            C(cnt) = K-sum1/sqrt(pi);
    end
    
    if length(find(C==R*K))==1
        fprintf("limit=%d\n",EbNo_dB(C==R*K));
        limits(R_cnt) = EbNo_dB(C==R*K);
        break;
    end
    
    tmp = find(C>R*K);
    g_index = tmp(1);
    l_index = tmp(1)-1;
    if(C(g_index)-C(l_index)<=1e-6)
        fprintf("limit=%d\n",EbNo_dB(l_index));
        limits(R_cnt) = EbNo_dB(g_index);
        break;
    end
    ebno_start = EbNo_dB(l_index);
    ebno_end = EbNo_dB(g_index);
    EbNo_dB = linspace(ebno_start,ebno_end,ebno_num);
    EbNo_Lin = 10.^(0.1*EbNo_dB);
    C = zeros(size(EbNo_dB));
end
end
toc
% figure;
% plot(limits,R_vec,'LineWidth',1);
% grid on;