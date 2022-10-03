% clear;

%% set simulation parameters
M = 256;% modulation order
K = log2(M);% bits per 2D
conste.dim = 2;% constellation dimension
conste.labels = 0:M-1;% decimal labels
conste.symbols = qammod(0:M-1, M, 'gray', 'UnitAveragePower', true);%getConste(M);% symbols
% conste.labels = qamdemod(conste.symbols, M,'gray', 'UnitAveragePower', true);
% p = 2.^[-6 -6 -5 -4 -4 -4 -3 -3 -3 -3 -4 -4 -4 -5 -6 -6];
p = ones(1,M)/M;
conste.p = p;
type = 'b';% CM('c') or BICM('b')

R = 1/2;% code rate
rho = R*K;% information bits per 2D (or spectra efficiency)


%% calculate limit
snr_start = 10;
snr_end = 14;
snr_num = 10;
SNR = linspace(snr_start,snr_end,snr_num);
while(1) 
    
    capacity = getCapacity_update(conste,SNR,type);

    tmp = find(capacity>rho);
    g_index = tmp(1);
    l_index = tmp(1)-1;
    if(capacity(g_index)-capacity(l_index)<=1e-2)
        limit = SNR(g_index);
        break;
    end
    
    snr_start = SNR(l_index);
    snr_end = SNR(g_index);
    fprintf("SNR_start: %d\n", snr_start)
    fprintf("SNR_end: %d\n", snr_end)
    SNR = linspace(snr_start,snr_end,snr_num);
    
end

% type = 'b';
% snr_start = 10;
% snr_end = 15;
% snr_num = 10;
% SNR = linspace(snr_start,snr_end,snr_num);
% while(1) 
%     
%     capacity = getCapacity_update(conste,SNR,type);
% 
%     tmp = find(capacity>rho);
%     g_index = tmp(1);
%     l_index = tmp(1)-1;
%     if(capacity(g_index)-capacity(l_index)<=1e-6)
%         limit_bicm = SNR(g_index);
%         break;
%     end
%     
%     snr_start = SNR(l_index);
%     snr_end = SNR(g_index);
%     SNR = linspace(snr_start,snr_end,snr_num);
%     
% end
%% output
% fprintf("Shannon limit: %d\n", 10*log10(2^rho-1))
% fprintf("The capacity limit is: \n");
% fprintf("CM-AMI limit: %d\n", limit);
fprintf("BICM limit: %d\n", limit_bicm);
% fprintf("    EsNo       %d\n", limit+10*log10(0.5*conste.dim));
% fprintf("    EbNo       %d\n", limit+10*log10(0.5*conste.dim/rho));