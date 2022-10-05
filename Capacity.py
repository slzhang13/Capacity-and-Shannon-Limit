import numpy as np
import torch
from tqdm import tqdm
from scipy import interpolate


class SISO_Capacity:
    def __init__(self, args_dict):

        self.symbols = args_dict["symbols"]
        self.chan_type = args_dict["chan_type"]

        self.M = self.symbols.size
        self.syms_dim = 2 if np.imag(self.symbols).any() else 1

        self.device = torch.device("cuda:0")

        self.syms_dist = torch.tensor(
            np.kron(self.symbols, np.ones((1, self.M)))
            - np.kron(np.ones((1, self.M)), self.symbols),
            dtype=torch.complex64,
            device=self.device,
        )

    def get_capacity(self, SNR_dB_list, SNR_type, sim_num):

        capacity = [0] * len(SNR_dB_list)
        SNR_bias = 0
        if SNR_type == "EsN0" and self.syms_dim == 1:
            SNR_bias = 10 * np.log10(2)
        if SNR_type == "EbN0":
            if self.syms_dim == 1:
                SNR_bias = 10 * np.log10(2 * np.log2(self.M))
            else:
                SNR_bias = 10 * np.log10(1 * np.log2(self.M))

        for i, SNR_dB in enumerate(SNR_dB_list):
            SNR_dB += SNR_bias
            capacity[i] = self.get_capacity_SNR(SNR_dB, sim_num)
        return capacity

    def get_capacity_SNR(self, SNR_dB, sim_num):

        SNR = 10 ** (0.1 * SNR_dB)
        sigma2 = 1 / SNR

        if self.syms_dim == 1:
            awgn = np.sqrt(sigma2) * torch.randn(
                sim_num, 1, dtype=torch.float32, device=self.device
            )
        else:
            awgn = np.sqrt(sigma2 / 2) * (
                torch.randn(sim_num, 1, dtype=torch.float32, device=self.device)
                + 1j * torch.randn(sim_num, 1, dtype=torch.float32, device=self.device)
            )

        if self.chan_type == "awgn":
            h = 1
        else:
            h = np.sqrt(1 / 2) * (
                torch.randn(sim_num, 1, dtype=torch.float32, device=self.device)
                + 1j * torch.randn(sim_num, 1, dtype=torch.float32, device=self.device)
            )

        if self.syms_dim == 1:
            metric = torch.exp(
                -abs(awgn + h * self.syms_dist).reshape(sim_num, self.M, self.M) ** 2
                / (2 * sigma2)
            )
        else:
            metric = torch.exp(
                -abs(awgn + h * self.syms_dist).reshape(sim_num, self.M, self.M) ** 2
                / sigma2
            )

        tmp = torch.log2(metric.sum(axis=2))
        tmp = torch.mean(tmp).cpu().numpy()

        return np.log2(self.M) - tmp

