# 1. 2-D constellation case

Suppose AWGN with $N_o=1$ for derivation simplicity.
$$
\begin{split}
C&=I(X;Y)\\
&=\sum_{x\in X}\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}p(x,y)\log\frac{p(y\mid x)}{p(y)}dy_{I}dy_{Q}\\
&=\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}p(y\mid x)\log\frac{p(y\mid x)}{\sum_{x^\prime\in X}p(y\mid x^\prime)p(x^\prime)}dy_{I}dy_{Q}\\
&=\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}\frac{1}{\pi}e^{-|y-x|^2}\log\frac{e^{-|y-x|^2}}{\sum_{x^\prime\in X}e^{-|y-x^\prime|^2}p(x^\prime)}dy_{I}dy_{Q}\\
&=-\frac{1}{\pi}\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}e^{-|y-x|^2}\log\frac{\sum_{x^\prime\in X}e^{-|y-x^\prime|^2}p(x^\prime)}{e^{-|y-x|^2}}dy_{I}dy_{Q}\\
&=-\frac{1}{\pi}\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}e^{-t_I^2}e^{-t_Q^2}\log{\sum_{x^\prime\in X}e^{|t|^2-|x-x^\prime+t|^2}p(x^\prime)}dt_{I}dt_{Q}\\
&=-\frac{1}{\pi}\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}e^{-t_I^2}e^{-t_Q^2}\log{\sum_{x^\prime\in X}e^{|t_I|^2+|t_Q|^2-|x-x^\prime+t_I+i\cdot t_Q|^2}p(x^\prime)}dt_{I}dt_{Q}\\
\end{split}\tag{1-1}
$$
When the constellation symbols are taken with equal probability, the above formula reduces to
$$
\begin{split}
C&=I(X;Y)\\
&=m-\frac{1}{M\pi}\sum_{x\in X}\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}e^{-t_I^2}e^{-t_Q^2}\log{\sum_{x^\prime\in X}e^{|t_I|^2+|t_Q|^2-|x-x^\prime+t_I+i\cdot t_Q|^2}}dt_{I}dt_{Q}\\
\end{split}\tag{1-2}
$$
where $m=\log_2M$ represents the number of coded bits per constellation symbol.

> Notice the calculation of average symbol energy when the constellation symbols are not taken with equal probability.
>
> $\bar{E_s}=\sum_{s\in \chi}p(s)|s|^2$

------



# 2. 1-D constellation case

Suppose AWGN with $N_o=1$ for derivation simplicity.
$$
\begin{split}
C&=I(X;Y)\\
&=\sum_{x\in X}\int_{-\infty}^{+\infty}p(x,y)\log\frac{p(y\mid x)}{p(y)}dy\\
&=\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}p(y\mid x)\log\frac{p(y\mid x)}{\sum_{x^\prime\in X}p(y\mid x^\prime)p(x^\prime)}dy\\
&=\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}\frac{1}{\sqrt{\pi}}e^{-|y-x|^2}\log\frac{e^{-|y-x|^2}}{\sum_{x^\prime\in X}e^{-|y-x^\prime|^2}p(x^\prime)}dy\\
&=-\frac{1}{\sqrt{\pi}}\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}e^{-|y-x|^2}\log\frac{\sum_{x^\prime\in X}e^{-|y-x^\prime|^2}p(x^\prime)}{e^{-|y-x|^2}}dy\\
&=-\frac{1}{\sqrt{\pi}}\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}e^{-|t|^2}\log{\sum_{x^\prime\in X}e^{|t|^2-|x-x^\prime+t|^2}p(x^\prime)}dt\\
&=-\frac{1}{\sqrt{\pi}}\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}e^{-t^2}\log{\sum_{x^\prime\in X}e^{t^2-(x-x^\prime+t)^2}p(x^\prime)}dt\\
\end{split}\tag{2-1}
$$
When the constellation symbols are taken with equal probability, the above formula reduces to
$$
\begin{split}
C&=I(X;Y)\\
&=m-\frac{1}{M\sqrt{\pi}}\sum_{x\in X}\int_{-\infty}^{+\infty}e^{-t^2}\log{\sum_{x^\prime\in X}e^{t^2-(x-x^\prime+t)^2}}dt\\
\end{split}\tag{2-2}
$$
where $m=\log_2M$ represents the number of coded bits per constellation symbol.

> Notice the calculation of average symbol energy when the constellation symbols are not taken with equal probability.
>
> $\bar{E_s}=\sum_{s\in \chi}p(s)|s|^2$

------



# 3. 1-D case is a special 2-D case 

$$
\begin{split}
C&=I(X;Y)\\
&=-\frac{1}{\pi}\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}e^{-t_I^2}e^{-t_Q^2}\log{\sum_{x^\prime\in X}e^{|t_I|^2+|t_Q|^2-|x-x^\prime+t_I+i\cdot t_Q|^2}p(x^\prime)}dt_{I}dt_{Q}\\
\end{split}\tag{3-1}
$$

When the constellation is 1-D, i.e. $x$ and $x^\prime$ is real, equation (3-1) can be rearranged.

$$
\begin{split}
C&=I(X;Y)\\
&=-\frac{1}{\pi}\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}e^{-t_I^2}e^{-t_Q^2}\log{\sum_{x^\prime\in X}e^{|t_I|^2+|t_Q|^2-|x-x^\prime+t_I|^2+|i\cdot t_Q|^2}p(x^\prime)}dt_{I}dt_{Q}\\

&=-\frac{1}{\pi}\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}e^{-t_I^2}e^{-t_Q^2}\log{\sum_{x^\prime\in X}e^{|t_I|^2-|x-x^\prime+t_I|^2}p(x^\prime)}dt_{I}dt_{Q}\\

&=-\frac{1}{\sqrt{\pi}}\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}e^{-t_I^2}\log{\sum_{x^\prime\in X}e^{|t_I|^2-|x-x^\prime+t_I|^2}p(x^\prime)}dt_{I}\int_{-\infty}^{+\infty}\frac{1}{\sqrt{\pi}}e^{-t_Q^2}dt_{Q}\\

&=-\frac{1}{\sqrt{\pi}}\sum_{x\in X}p(x)\int_{-\infty}^{+\infty}e^{-t_I^2}\log{\sum_{x^\prime\in X}e^{|t_I|^2-|x-x^\prime+t_I|^2}p(x^\prime)}dt_{I}
\end{split}\tag{3-2}
$$

Equation (3-2) is just an equivalent form of equation (2-1).

------



# 4. Relationship between several form of signal to noise ratio

$$
SNR=\frac{\sigma^2_s}{\sigma_n^2}=\begin{cases}
\frac{E_s}{N_o}\cdot2 &\quad\quad \text{1-D constellation}\\
\\
\\
\frac{E_s}{N_o} &\quad\quad \text{2-D constellation}\\
\end{cases}\tag{4-1}
$$


$$
\frac{E_s}{N_o}=\frac{E_b*\rho}{N_o}
\tag{4-2}
$$
where $\rho$ is the data rate in bits per symbol (or channel use) and
$$
\rho=m\cdot R
\tag{4-3}
$$
where $m$ is the number of coded bits per constellation symbol and $R$ is the code rate of FEC.