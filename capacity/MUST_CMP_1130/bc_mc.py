import numpy as np
import torch

def de2bi(nums, K):
    
    if(nums.size>1):
        nums = nums.squeeze() #将nums尺寸统一为(M, )
    
    M = len(nums)

    ret = []
    for i in range(M): #遍历所有要转换的十进制数
        num = nums[i]
        tmp = []
        for i in range(K-1, -1, -1): #逐二进制位遍历
            tmp.append((num >> i) & 1)
        ret.append(tmp)

    ret = np.array(ret)
    
    if ret.shape[0]==1:
        ret = ret.squeeze()
        
    return ret

def awgn(symbols, labels, Kb=0):

    
    N_samples = int(1e4)
    
    M = symbols.size
    K = int(np.log2(M))
    if Kb==0:
        Kb=K
    if labels.shape[1] == 1: #若传入的标识为十进制，先转换为二进制
        labels = de2bi(labels, K)
        
    device = torch.device('cuda:0')
    
    symbols = torch.tensor(symbols, dtype=torch.complex64, device=device)
    
    n = np.sqrt(1/2)*(torch.randn(N_samples, 1, M, dtype=torch.float32, device=device) +
      1j*torch.randn(N_samples, 1, M, dtype=torch.float32, device=device))
        
    metric = torch.exp(-abs(n+symbols.T-symbols)**2)
    
    bit_capacity = np.zeros(Kb)

    for k in range(Kb):
        
        mask0 = np.array(np.ix_(labels[:, k]==0))
        mask1 = np.array(np.ix_(labels[:, k]==1))
    
        p00 = torch.sum(metric[:, mask0.T, mask0], 1)
        p10 = torch.sum(metric[:, mask1.T, mask0], 1)
        p01 = torch.sum(metric[:, mask0.T, mask1], 1)
        p11 = torch.sum(metric[:, mask1.T, mask1], 1)
        
        tmp0 = torch.mean(torch.log2(1+p10/p00))
        tmp1 = torch.mean(torch.log2(1+p01/p11))
        bit_capacity[k] = 1-0.5*(tmp0+tmp1)

    torch.cuda.empty_cache()
            
    return bit_capacity


def awgn_legacy(symbols1, labels1, P2, Kb=0):

    
    N_samples = int(1e4)
    
    M1 = symbols1.size
    K1 = int(np.log2(M1))
    if Kb==0:
        Kb=K1
    if labels1.shape[1] == 1: #若传入的标识为十进制，先转换为二进制
        labels1 = de2bi(labels1, K1)
        
    device = torch.device('cuda:0')
    
    symbols1 = torch.tensor(symbols1, dtype=torch.complex64, device=device)
    
    n = np.sqrt((1+P2)/2)*(torch.randn(N_samples, 1, M1, dtype=torch.float32, device=device) +
      1j*torch.randn(N_samples, 1, M1, dtype=torch.float32, device=device))
        
    metric = torch.exp(-abs(n+symbols1.T-symbols1)**2/(1+P2))
    
    bit_capacity = np.zeros(Kb)

    for k in range(Kb):
        
        mask10 = np.array(np.ix_(labels1[:, k]==0))
        mask11 = np.array(np.ix_(labels1[:, k]==1))
    
        p00 = torch.sum(metric[:, mask10.T, mask10], 1)
        p10 = torch.sum(metric[:, mask11.T, mask10], 1)
        p01 = torch.sum(metric[:, mask10.T, mask11], 1)
        p11 = torch.sum(metric[:, mask11.T, mask11], 1)
        
        tmp0 = torch.mean(torch.log2(1+p10/p00))
        tmp1 = torch.mean(torch.log2(1+p01/p11))
        bit_capacity[k] = 1-0.5*(tmp0+tmp1)

    torch.cuda.empty_cache()
            
    return bit_capacity

def rayleigh(symbols, labels, Kb=0):

    
    N_samples = int(1e4)
    
    M = symbols.size
    K = int(np.log2(M))
    if Kb==0:
        Kb=K   
    if labels.shape[1] == 1: #若传入的标识为十进制，先转换为二进制
        labels = de2bi(labels, K)
        
    device = torch.device('cuda:0')
    
    symbols = torch.tensor(symbols, dtype=torch.complex64, device=device)
    
    n = np.sqrt(1/2)*(torch.randn(N_samples, 1, M, dtype=torch.float32, device=device) +
      1j*torch.randn(N_samples, 1, M, dtype=torch.float32, device=device))
    
    h = np.sqrt(1/2)*(torch.randn(N_samples, 1, M, dtype=torch.float32, device=device) +
      1j*torch.randn(N_samples, 1, M, dtype=torch.float32, device=device))
        
    metric = torch.exp(-abs(n+h*(symbols.T-symbols))**2)
    
    bit_capacity = np.zeros(Kb)

    for k in range(Kb):
        
        mask0 = np.array(np.ix_(labels[:, k]==0))
        mask1 = np.array(np.ix_(labels[:, k]==1))
    
        p00 = torch.sum(metric[:, mask0.T, mask0], 1)
        p10 = torch.sum(metric[:, mask1.T, mask0], 1)
        p01 = torch.sum(metric[:, mask0.T, mask1], 1)
        p11 = torch.sum(metric[:, mask1.T, mask1], 1)
        
        tmp0 = torch.mean(torch.log2(1+p10/p00))
        tmp1 = torch.mean(torch.log2(1+p01/p11))
        bit_capacity[k] = 1-0.5*(tmp0+tmp1)

    torch.cuda.empty_cache()
            
    return bit_capacity


def rayleigh_legacy(symbols1, labels1, P2, Kb=0):

    
    N_samples = int(1e4)
    
    M1 = symbols1.size
    K1 = int(np.log2(M1))
    if Kb==0:
        Kb=K1    
    if labels1.shape[1] == 1: #若传入的标识为十进制，先转换为二进制
        labels1 = de2bi(labels1, K1)
        
    device = torch.device('cuda:0')

    symbols1 = torch.tensor(symbols1, dtype=torch.complex64, device=device)
    
    h = np.sqrt(1/2)*(torch.randn(N_samples, 1, M1, dtype=torch.float32, device=device) +
      1j*torch.randn(N_samples, 1, M1, dtype=torch.float32, device=device))
    
    n = torch.sqrt((1+abs(h)**2*P2)/2)*(torch.randn(N_samples, 1, M1, dtype=torch.float32, device=device) +
      1j*torch.randn(N_samples, 1, M1, dtype=torch.float32, device=device))
    
    metric = torch.exp(-abs(n+h*(symbols1.T-symbols1))**2/(1+abs(h)**2*P2))
    
    bit_capacity = np.zeros(Kb)

    for k in range(Kb):
        
        mask10 = np.array(np.ix_(labels1[:, k]==0))
        mask11 = np.array(np.ix_(labels1[:, k]==1))
    
        p00 = torch.sum(metric[:, mask10.T, mask10], 1)
        p10 = torch.sum(metric[:, mask11.T, mask10], 1)
        p01 = torch.sum(metric[:, mask10.T, mask11], 1)
        p11 = torch.sum(metric[:, mask11.T, mask11], 1)
        
        tmp0 = torch.mean(torch.log2(1+p10/p00))
        tmp1 = torch.mean(torch.log2(1+p01/p11))
        bit_capacity[k] = 1-0.5*(tmp0+tmp1)

    torch.cuda.empty_cache()
            
    return bit_capacity