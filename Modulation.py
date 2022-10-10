import numpy as np

round_digit = 6


def de2bi(dec_nums, K=0):

    # modulation order in current context
    M = dec_nums.size
    if K == 0:
        K = int(np.log2(M))

    ret = np.zeros((M, K), dtype=np.int32)

    # for each decimal label
    for m in range(M):
        # for each bit
        for k in range(K):  # left-msb
            ret[m, k] = int((dec_nums[0, m] >> (K - 1 - k)) & 1)

    return ret


def grayCoded(n):
    # return Gray-coded integer numbers 0~n-1
    ret = np.zeros((1, n), dtype=np.int32)

    for i in range(n):
        ret[0, i] = i ^ (i >> 1)

    return ret


def normalizeConste(symbols):

    M = symbols.size

    normalized_symbols = symbols / np.linalg.norm(symbols) * np.sqrt(M)

    return normalized_symbols


def PSK(M, labeling="gray"):
    # start from point (1, 0), anti-clockwise

    symbols = np.exp(1j * 2 * np.pi / M * np.arange(M)).reshape(1, M)
    labels = grayCoded(M)

    if labeling == "gray":
        pass
    elif labeling == "bin":
        labels = np.arange(M).reshape(1, M)
    else:
        print("Unsupported labeling name! Gray labeling is used instead.")
        print("But you can customize the labeling with a vector outside this function.")

    symbols = symbols.round(round_digit)  # 1-D/2-D justification！！！

    return symbols, labels


def QAM(M, labeling="gray"):
    # rightwards and downwards

    if isinstance(M, int):
        MI = MQ = int(np.sqrt(M))
    elif isinstance(M, tuple):
        MI, MQ = M[0], M[1]
        M = MI * MQ
    else:
        pass

    I_symbols = np.array([np.arange(-MI + 1, MI + 1, 2)])
    Q_symbols = -1j * np.array([np.arange(-MQ + 1, MQ + 1, 2)]).T

    symbols = I_symbols + Q_symbols

    I_labels = grayCoded(MI)  # 维度为(1, MI)
    Q_labels = grayCoded(MQ).T * MI  # 维度为(MI, 1)

    labels = I_labels + Q_labels  # 维度为(MI, MI)

    symbols = symbols.reshape(1, M)
    labels = labels.reshape(1, M)

    symbols = normalizeConste(symbols)  # 归一化能量

    if labeling == "gray":
        pass
    elif labeling == "bin":
        labels = np.arange(M).reshape(1, M)
    else:
        print("Unsupported labeling name! Gray labeling is used instead.")
        print("But you can customize the labeling with a vector outside this function.")

    symbols = symbols.round(round_digit)  # 1-D/2-D justification！！！

    return symbols, labels


def PAM(M, labeling="gray"):
    # rightwards

    symbols = np.array([np.arange(-M + 1, M + 1, 2)])
    labels = grayCoded(M)

    symbols = symbols.reshape(1, M)
    labels = labels.reshape(1, M)

    symbols = normalizeConste(symbols)  # 归一化能量

    if labeling == "gray":
        pass
    elif labeling == "bin":
        labels = np.arange(M).reshape(1, M)
    else:
        print("Unsupported labeling name! Gray labeling is used instead.")
        print("But you can customize the labeling with a vector outside this function.")

    symbols = symbols.round(round_digit)  # 1-D/2-D justification！！！

    return symbols, labels

