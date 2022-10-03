M = 16;
conste.dim = 2;
conste.symbols = qammod(0:M-1, M,'UnitAveragePower', true);
conste.labels = 0:M-1;
SNR = 10;
type = 'c';
get_capacity(conste,SNR,type)
cal_capacity(conste.symbols,conste.dim,SNR)