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
