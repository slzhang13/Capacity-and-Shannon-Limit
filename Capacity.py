## import necessary packages
import numpy as np
import torch
from tqdm import tqdm
from Modulation import de2bi


## Gaussian Inputs
def Gaussian_Capacity(chan_type, SNR_range, sim_num=10000):

    if chan_type == "awgn":

        capacity = np.log2(1 + 10 ** (0.1 * SNR_range))

    elif chan_type == "rayleigh":

        capacity = np.zeros(SNR_range.shape)

        # GPU device
        device = torch.device("cuda:0")

        h = np.sqrt(1 / 2) * (
            torch.randn(1, sim_num, dtype=torch.float32, device=device)
            + 1j * torch.randn(1, sim_num, dtype=torch.float32, device=device)
        )
        for i, SNR_dB in tqdm(enumerate(SNR_range)):
            SNR = 10 ** (0.1 * SNR_dB)
            metric = torch.log2(1 + abs(h) ** 2 * SNR)
            capacity[i] = metric.mean().cpu().numpy()

    else:
        print("Unsupported channel type!")
        exit()

    return capacity


## Finite Alphabet Inputs (CM)
def CM_Capacity(conste_symbols, chan_type, SNR_range, ssd=(False, 0), sim_num=10000):

    # constellation rotation in signal space diversity
    if ssd[0]:
        conste_symbols *= np.exp(1j * ssd[1] / 180 * np.pi)

    # modulation order
    M = conste_symbols.size

    # dimension of constellation, 2 for I/Q modulation, and 1 for BPSK/PAM
    sym_dim = 2 if np.imag(conste_symbols).any() else 1

    # GPU device
    device = torch.device("cuda:0")

    # euclidean distances between constellation points
    sym_dists = torch.tensor(
        np.kron(conste_symbols, np.ones((1, M)))
        - np.kron(np.ones((1, M)), conste_symbols),
        dtype=torch.complex64,
        device=device,
    )

    # output capacity results
    capacity = np.zeros(SNR_range.shape)

    # awgn realizations (0.5 variance for each noise dimmension)
    if sym_dim == 1:
        awgn = np.sqrt(0.5) * torch.randn(
            sim_num, 1, dtype=torch.float32, device=device
        )
    else:
        awgn = np.sqrt(0.5) * (
            torch.randn(sim_num, 1, dtype=torch.float32, device=device)
            + 1j * torch.randn(sim_num, 1, dtype=torch.float32, device=device)
        )

    # channel realizations
    if chan_type == "awgn":
        hI = hQ = 1
    elif chan_type == "rayleigh":
        hI = np.sqrt(1 / 2) * (
            torch.randn(sim_num, 1, dtype=torch.float32, device=device)
            + 1j * torch.randn(sim_num, 1, dtype=torch.float32, device=device)
        )
        if ssd[0]:
            hQ = np.sqrt(1 / 2) * (
                torch.randn(sim_num, 1, dtype=torch.float32, device=device)
                + 1j * torch.randn(sim_num, 1, dtype=torch.float32, device=device)
            )
        else:
            hQ = hI
    else:
        print("Unsupported channel type!")
        exit()

    # main loop
    for i, SNR_dB in tqdm(enumerate(SNR_range)):

        SNR = 10 ** (0.1 * SNR_dB)

        # signal power
        Ps = SNR * sym_dim * 0.5

        sym_dists_I = np.sqrt(Ps) * sym_dists.real
        sym_dists_Q = np.sqrt(Ps) * sym_dists.imag

        metric = torch.exp(
            -abs(awgn + abs(hI) * sym_dists_I + abs(hQ) * sym_dists_Q).reshape(
                sim_num, M, M
            )
            ** 2
        )

        tmp = torch.log2(metric.sum(axis=2))
        tmp = torch.mean(tmp).cpu().numpy()

        capacity[i] = np.log2(M) - tmp - 0.5 * sym_dim * np.log2(np.exp(1))

    return capacity


## Finite Alphabet Inputs (BICM)
def BICM_Capacity(
    conste_symbols, conste_labels, chan_type, SNR_range, ssd=(False, 0), sim_num=10000
):
    # constellation rotation in signal space diversity
    if ssd[0]:
        conste_symbols *= np.exp(1j * ssd[1] / 180 * np.pi)

    conste_labels = de2bi(conste_labels)

    # modulation order
    M = conste_symbols.size

    # bits per constellation point
    K = int(np.log2(M))

    # dimension of constellation, 2 for I/Q modulation, and 1 for BPSK/PAM
    sym_dim = 2 if np.imag(conste_symbols).any() else 1

    # GPU device
    device = torch.device("cuda:0")

    # euclidean distances between constellation points
    sym_dists = torch.tensor(
        np.kron(conste_symbols, np.ones((1, M)))
        - np.kron(np.ones((1, M)), conste_symbols),
        dtype=torch.complex64,
        device=device,
    )

    # output capacity results
    bit_capacity = np.zeros((SNR_range.size, K))

    # awgn realizations (0.5 variance for each noise dimmension)
    if sym_dim == 1:
        awgn = np.sqrt(0.5) * torch.randn(
            sim_num, 1, dtype=torch.float32, device=device
        )
    else:
        awgn = np.sqrt(0.5) * (
            torch.randn(sim_num, 1, dtype=torch.float32, device=device)
            + 1j * torch.randn(sim_num, 1, dtype=torch.float32, device=device)
        )

    if chan_type == "awgn":
        hI = hQ = 1
    elif chan_type == "rayleigh":
        hI = np.sqrt(1 / 2) * (
            torch.randn(sim_num, 1, dtype=torch.float32, device=device)
            + 1j * torch.randn(sim_num, 1, dtype=torch.float32, device=device)
        )
        if ssd[0]:
            hQ = np.sqrt(1 / 2) * (
                torch.randn(sim_num, 1, dtype=torch.float32, device=device)
                + 1j * torch.randn(sim_num, 1, dtype=torch.float32, device=device)
            )
        else:
            hQ = hI
    else:
        print("Unsupported channel type!")
        exit()

    # main loop
    for i, SNR_dB in tqdm(enumerate(SNR_range)):

        SNR = 10 ** (0.1 * SNR_dB)

        # signal power
        Ps = SNR * sym_dim * 0.5

        sym_dists_I = np.sqrt(Ps) * sym_dists.real
        sym_dists_Q = np.sqrt(Ps) * sym_dists.imag

        metric = torch.exp(
            -abs(awgn + abs(hI) * sym_dists_I + abs(hQ) * sym_dists_Q).reshape(
                sim_num, M, M
            )
            ** 2
        )

        for k in range(K):

            mask0 = np.array(np.ix_(conste_labels[:, k] == 0))
            mask1 = np.array(np.ix_(conste_labels[:, k] == 1))

            p00 = torch.sum(metric[:, mask0.T, mask0], 2)
            p01 = torch.sum(metric[:, mask0.T, mask1], 2)
            p10 = torch.sum(metric[:, mask1.T, mask0], 2)
            p11 = torch.sum(metric[:, mask1.T, mask1], 2)

            tmp0 = torch.mean(torch.log2(1 + p01 / p00)).cpu().numpy()
            tmp1 = torch.mean(torch.log2(1 + p10 / p11)).cpu().numpy()
            bit_capacity[i, k] = 1 - 0.5 * (tmp0 + tmp1)

    capacity = bit_capacity.sum(axis=1)

    return bit_capacity, capacity


def getCapacity(
    conste_symbols,
    conste_labels,
    chan_type,
    CM_type,
    SNR_range,
    ssd=(False, 0),
    sim_num=10000,
):
    if CM_type == "gaussian":
        return Gaussian_Capacity(chan_type, SNR_range, sim_num)
    elif CM_type == "CM":
        return CM_Capacity(conste_symbols, chan_type, SNR_range, ssd, sim_num)
    else:
        return BICM_Capacity(
            conste_symbols, conste_labels, chan_type, SNR_range, ssd, sim_num,
        )

