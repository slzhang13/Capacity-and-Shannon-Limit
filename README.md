# Capacity-and-Shannon-Limit
This repository provides MATLAB codes to calculate the capacity and shannon limit for a specific coded modulation scheme.

#### Relationships between different SNR definitions

$$
\text{SNR}=\frac{\sigma_s^2}{\sigma_n^2}=\frac{E_s}{N_0}\times \{1, 2\}
$$
>2 for 1-D constellations (BPSK, PAM)
1 for 2-D constellations


$$
\frac{E_s}{N_0}=\frac{E_b(K\times R)}{N_0}
$$
>K denotes number of (coded) bits per constellation symbol
R denotes the code rate

$$
\text{SNR}=\frac{E_b}{N_0}\times \eta
$$

>$\eta$ can be interpreted as the number of information bits per 2-D signal

#### Capacity VS. $E_b/N_0$
+ capacity vs. SNR
+ using the fact that $E_b/N_0$ is a function of SNR and capacity
> For Gaussian inputs, there is a closed form expression:
$$\frac{E_b}{N_0}=\frac{2^\eta-1}{\eta}$$

#### capacity for Gaussian inputs:
+ At low SNR, $C_{\text{Rayleigh}}\approx C_\text{AWGN}$
+ At high SNR, $C_{\text{Rayleigh}}\approx C_\text{AWGN}+0.83$
> TSE Book P235
0.83 comes from the mean of Log-Rayleigh Distribution
