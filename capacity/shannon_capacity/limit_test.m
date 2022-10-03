clear;
M = 2;
conste_norm = exp(1j*(0:M-1)*2*pi/M);
% conste_tmp = [-3 -1 1 3]+1j*[3; 1; -1; -3];
% conste_tmp = conste_tmp(:);
% conste_norm = conste_tmp./sqrt(mean(abs(conste_tmp).^2));
K = log2(M);
R_vec = 1/2;
limits = zeros(size(R_vec));
for R_cnt = 1:length(R_vec)
    R = R_vec(R_cnt);
    ebno_start = -20;
    ebno_end = 40;
    ebno_num = 100;
    ebno_dB = linspace(ebno_start,ebno_end,ebno_num);
    
    iter = 1;
    while(1)
        fprintf("%d\n",iter);
        iter = iter+1;
        C = cal_capacity(conste_norm, R, ebno_dB);
        
        if length(find(C==R*K))==1
            fprintf("limit=%d\n",ebno_dB(C==R*K));
            limits(R_cnt) = ebno_dB(C==R*K);
            break;
        end

        tmp = find(C>R*K);
        g_index = tmp(1);
        l_index = tmp(1)-1;
        if(C(g_index)-C(l_index)<=1e-6)
            fprintf("limit=%d\n",ebno_dB(l_index));
            limits(R_cnt) = ebno_dB(l_index);
            break;
        end
        ebno_start = ebno_dB(l_index);
        ebno_end = ebno_dB(g_index);
        ebno_dB = linspace(ebno_start,ebno_end,ebno_num);
        C = zeros(size(ebno_dB));
    end
end