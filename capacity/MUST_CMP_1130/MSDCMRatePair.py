from typing_extensions import runtime
import numpy as np
import torch

from bc_gh import de2bi


def getMSDCMRatePair(symbols1, labels1, symbols2, labels2, SNR1_lin, SNR2_lin, alpha):

    M1 = symbols1.size
    M2 = symbols2.size
    K1 = int(np.log2(M1))
    K2 = int(np.log2(M2))

    K = K1+K2
    M = M1*M2

    device = torch.device('cuda:0')

    xx, ww = GH(n=16)
    xx = xx.reshape((xx.size, 1, 1))
    xx = torch.tensor(xx, dtype=torch.complex64, device=device)
    ww = torch.tensor(ww, dtype=torch.float32, device=device)

    # 基本层容量推导
    symbols = torch.sqrt(SNR1_lin)*torch.tensor(symbols1*np.sqrt(1-alpha)+symbols2.T*np.sqrt(alpha),
                           dtype=torch.complex64, device=device)
    symbols = np.reshape(symbols, (M, 1), 'F')

    labels = labels1*M2+labels2.T
    labels = np.reshape(labels, (M, 1), 'F')

    labels = de2bi(labels, K)
    metric = torch.exp(-abs(xx-symbols+symbols.T)**2)

    Rbc = 0

    for k in range(K1):

        mask0 = np.array(np.ix_(labels[:, k] == 0))
        mask1 = np.array(np.ix_(labels[:, k] == 1))

        p00 = torch.sum(metric[:, mask0.T, mask0], 1)
        p10 = torch.sum(metric[:, mask1.T, mask0], 1)
        p01 = torch.sum(metric[:, mask0.T, mask1], 1)
        p11 = torch.sum(metric[:, mask1.T, mask1], 1)

        Rbc += 1-torch.sum(ww*torch.sum(torch.log2(
            1+p10/p00)+torch.log2(1+p01/p11), axis=1))/(M*np.pi)

    del metric

    # 增强层
    s2 = torch.sqrt(SNR2_lin)*torch.tensor(symbols2*np.sqrt(alpha),
                      dtype=torch.complex64, device=device)
    labels2 = de2bi(labels2, K2)

    metric = torch.exp(-abs(xx-s2+s2.T)**2)

    Ruc = 0

    for k in range(K2):

        mask0 = np.array(np.ix_(labels2[:, k] == 0))
        mask1 = np.array(np.ix_(labels2[:, k] == 1))

        p00 = torch.sum(metric[:, mask0.T, mask0], 1)
        p10 = torch.sum(metric[:, mask1.T, mask0], 1)
        p01 = torch.sum(metric[:, mask0.T, mask1], 1)
        p11 = torch.sum(metric[:, mask1.T, mask1], 1)

        Ruc += 1-torch.sum(ww*torch.sum(torch.log2(
            1+p10/p00)+torch.log2(1+p01/p11), axis=1))/(M2*np.pi)

    return Rbc, Ruc


def getMSDBICMRatePair(symbols1, symbols2, SNR1_lin, SNR2_lin, alpha):

    M1 = symbols1.size
    M2 = symbols2.size

    M = M1*M2

    device = torch.device('cuda:0')

    xx, ww = GH(n=16)
    xx = xx.reshape((xx.size, 1, 1))
    xx = torch.tensor(xx, dtype=torch.complex64, device=device)
    ww = torch.tensor(ww, dtype=torch.float32, device=device)

    # 基本层
    s1 = torch.tensor(symbols1*np.sqrt(1-alpha)),
                      dtype = torch.complex64, device = device)
    s2=torch.tensor(symbols2*np.sqrtalpha),
                      dtype=torch.complex64, device=device)

    symbols=s1+s2.T
    symbols=symbols.reshape((M, 1))

    metric_num=torch.exp(-abs(xx-symbols+symbols.T)**2)
    metric_den=torch.exp(-abs(xx-s2+s2.T)**2)

    tmp_num=torch.log2(torch.sum(metric_num, axis = 1))
    tmp_den=torch.log2(torch.sum(metric_den, axis = 1))

    tmp1=torch.sum(ww*torch.sum(tmp_num, axis = 1))/(M*np.pi)
    tmp2=torch.sum(ww*torch.sum(tmp_den, axis = 1))/(M2*np.pi)

    Rbc=np.log2(M1)-tmp1+tmp2

    # 增强层
    s2=torch.tensor(symbols2*np.sqrt(SNR2_lin*alpha),
                      dtype=torch.complex64, device=device)

    metric2=torch.exp(-abs(xx-s2+s2.T)**2+abs(xx)**2)

    tmp=torch.log2(torch.sum(metric2, axis = 1))

    Ruc=np.log2(M2)-torch.sum(ww*torch.sum(tmp, axis = 1))/(M2*np.pi)

    return Rbc, Ruc


def awgn(symbols, labels, Kb=0):

    M=symbols.size
    K=int(np.log2(M))

    if Kb == 0:
        Kb=K
        Krange=range(Kb)
    elif Kb > 0:
        Krange=range(Kb)
    else:
        Kb=-Kb
        Krange=range(K-Kb, K)

    if labels.shape[1] != K:
        labels=de2bi(labels, K)

    device=torch.device('cuda:0')

    xx, ww=GH(n=16)
    xx=xx.reshape((xx.size, 1, 1))
    xx=torch.tensor(xx, dtype=torch.complex64, device=device)
    ww=torch.tensor(ww, dtype=torch.float32, device=device)

    symbols=torch.tensor(symbols, dtype=torch.complex64, device=device)

    metric=torch.exp(-abs(xx-symbols+symbols.T)**2)

    bit_capacity=np.zeros(Kb)

    for i in range(len(Krange)):

        k=Krange[i]
        mask0=np.array(np.ix_(labels[:, k] == 0))
        mask1=np.array(np.ix_(labels[:, k] == 1))

        p00=torch.sum(metric[:, mask0.T, mask0], 1)
        p10=torch.sum(metric[:, mask1.T, mask0], 1)
        p01=torch.sum(metric[:, mask0.T, mask1], 1)
        p11=torch.sum(metric[:, mask1.T, mask1], 1)

        bit_capacity[i]=1-torch.sum(ww*torch.sum(torch.log2(
            1+p10/p00)+torch.log2(1+p01/p11), axis = 1))/(M*np.pi)

    return bit_capacity


def awgn_legacy(symbols, labels, P2, Kb=0):

    M=symbols.size
    K=int(np.log2(M))
    if Kb == 0:
        Kb=K
    if labels.shape[1] != K:
        labels=de2bi(labels, K)

    device=torch.device('cuda:0')

    xx, ww=GH(n=16)
    xx=xx.reshape((xx.size, 1, 1))
    xx=torch.tensor(xx, dtype=torch.complex64, device=device)
    ww=torch.tensor(ww, dtype=torch.float32, device=device)

    symbols=torch.tensor(symbols/np.sqrt(1+P2),
                           dtype=torch.complex64, device=device)

    metric=torch.exp(-abs(xx-symbols+symbols.T)**2)

    bit_capacity=np.zeros(Kb)

    for k in range(Kb):

        mask0=np.array(np.ix_(labels[:, k] == 0))
        mask1=np.array(np.ix_(labels[:, k] == 1))

        p00=torch.sum(metric[:, mask0.T, mask0], 1)
        p10=torch.sum(metric[:, mask1.T, mask0], 1)
        p01=torch.sum(metric[:, mask0.T, mask1], 1)
        p11=torch.sum(metric[:, mask1.T, mask1], 1)

        bit_capacity[k]=1-torch.sum(ww*torch.sum(torch.log2(
            1+p10/p00)+torch.log2(1+p01/p11), axis = 1))/(M*np.pi)

    return bit_capacity


def de2bi(nums, K):

    if(nums.size > 1):
        nums=nums.squeeze()  # 将nums尺寸统一为(M, )

    M=len(nums)

    ret=[]
    for i in range(M):  # 遍历所有要转换的十进制数
        num=nums[i]
        tmp=[]
        for i in range(K-1, -1, -1):  # 逐二进制位遍历
            tmp.append((num >> i) & 1)
        ret.append(tmp)

    ret=np.array(ret)

    if ret.shape[0] == 1:
        ret=ret.squeeze()

    return ret


def GH(n=16):

    p=np.zeros((n+1, n+1))
    p[0, n]=1
    p[1, n-1:]=np.array([1, 0])*2
    for k in range(2, n+1):
        p[k, n-k:]=2*np.append(p[k-1, n-k+1:], 0) - \
            2*(k-1)*np.insert(p[k-2, n-k+2:], 0, [0, 0])

    x=np.roots(p[n, :])
    w=np.zeros(n)
    for k in range(n):
        w[k]=2**(n-1)*np.math.factorial(n)*np.sqrt(np.pi) / \
            (n**2*np.polyval(p[n-1, :], x[k])**2)

    x=np.reshape(x, (1, n))
    w=np.reshape(w, (1, n))

    xx=x+1j*x.T
    ww=w*w.T

    xx=xx.flatten(order='F')
    ww=ww.flatten(order='F')

    xx=xx.reshape((1, len(xx)))
    ww=ww.reshape((1, len(ww)))

    return xx, ww
