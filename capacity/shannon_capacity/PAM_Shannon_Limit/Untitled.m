R=1e-10:1e-5:2;
ebno_limit = 10*log10((2.^(2*R)-1)./(2*R));
plot(ebno_limit,R)
grid on