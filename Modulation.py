import numpy as np

# all 1-d numpy arrays are row vectors


def de2bi(dec_nums, K):

    # modulation order in current context
    M = dec_nums.size

    ret = np.zeros((M, K))

    # for each decimal label
    for m in range(M):
        # for each bit
        for k in range(K):  # left-msb
            ret[m, k] = (dec_nums[m] >> (K - 1 - k)) & 1

    return ret


def grayCoded(n):
    # return Gray-coded integer numbers 0~n-1
    ret = np.zeros(n, dtype=np.int32)

    for i in range(n):
        ret[i] = i ^ (i >> 1)

    return ret.reshape(1, n)


def normalizeConste(symbols, sym_probs=0):

    M = symbols.size
    if sym_probs == 0:
        sym_probs = np.array([1 / M] * M).reshape(1, M)

    normalized_symbols = symbols / np.sqrt(abs(symbols) ** 2 @ sym_probs.T)

    return normalized_symbols


def PSK(M, labeling="gray"):

    symbols = np.exp(1j * 2 * np.pi / M * np.arange(M)).reshape(1, M)
    if labeling == "gray":
        labels = grayCoded(M).reshape(1, M)
    elif labeling == "bin":
        labels = np.arange(M).reshape(1, M)
    else:
        print("Unsupported labeling name! Gray labeling is used instead.")
        print("But you can customize the labeling with a vector outside this function.")
        labels = grayCoded(M).reshape(1, M)

    return symbols, labels
