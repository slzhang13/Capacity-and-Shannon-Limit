import numpy as np
import torch
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist
from scipy import interpolate


def de2bi(nums, K):

    if nums.size > 1:
        nums = nums.squeeze()  # 将nums尺寸统一为(M, )

    M = len(nums)

    ret = []
    for i in range(M):  # 遍历所有要转换的十进制数
        num = nums[i]
        tmp = []
        for i in range(K - 1, -1, -1):  # 逐二进制位遍历
            tmp.append((num >> i) & 1)
        ret.append(tmp)

    ret = np.array(ret)

    if ret.shape[0] == 1:
        ret = ret.squeeze()

    return ret


def grayCoded(n):

    ret = np.zeros((1, n), dtype=np.int32)

    for i in range(n):
        ret[0, i] = i ^ (i >> 1)

    return ret


def normalizeConste(symbols):

    M = len(symbols)
    normalized_symbols = symbols / np.linalg.norm(symbols) * np.sqrt(M)

    return normalized_symbols


def adjustBinLabel(symbols, labels):

    M = symbols.size
    K = int(np.log2(M))

    tmp1 = symbols.real > 0
    tmp2 = symbols.imag > 0
    tmp = np.zeros(M, dtype=bool)

    for i in range(M):
        tmp[i] = tmp1[i, 0] and tmp2[i, 0]

    labels_first_quad = labels[tmp, :]

    idx = []
    for k in range(K):
        if labels_first_quad[:, k].all() or (not labels_first_quad[:, k].any()):
            idx.append(k)

    labels[:, [0, idx[0]]] = labels[:, [idx[0], 0]]
    labels[:, [1, idx[1]]] = labels[:, [idx[1], 1]]

    return labels


def grayQAM(M, dec_label=False):

    K = int(np.log2(M))
    MI = int(np.sqrt(M))
    I_symbols = np.array([np.arange(-MI + 1, MI + 1, 2)])
    Q_symbols = -1j * I_symbols.T

    symbols = I_symbols + Q_symbols

    I_labels = grayCoded(MI)  # 维度为(1, MI)
    Q_labels = I_labels.T * MI  # 维度为(MI, 1)

    labels = I_labels + Q_labels  # 维度为(MI, MI)

    symbols = symbols.reshape((M, 1))
    labels = labels.reshape((M, 1))

    symbols = normalizeConste(symbols)  # 归一化能量

    labels = de2bi(labels, K)

    labels = adjustBinLabel(symbols, labels)  # 调整为两个高位比特决定象限

    if dec_label:  # 转换为二进制标识
        labels = labels @ (2 ** np.array(range(K - 1, -1, -1)))
        labels = labels.reshape((M, 1))

    return symbols, labels


def getAPSKRadius(m2):

    R_APSK = np.zeros((1, int(2 ** m2)))
    for i in range(R_APSK.size):
        R_APSK[0, i] = np.sqrt(-np.log(1 - (i + 1 - 0.5) * 2 ** (-m2)))

    return R_APSK


def grayAPSK(M_APSK, dec_label=False):

    M_APSK = M_APSK.squeeze()

    M = np.sum(M_APSK)  # 调制阶数
    K = int(np.log2(M))

    ring_num = M_APSK.size  # APSK环数
    m2 = int(np.log2(ring_num))

    # 参见参考文献中 Fig.1
    theta0 = np.pi / M_APSK[0]
    delta_theta = 2 * theta0

    R_APSK = getAPSKRadius(m2)
    thetas = np.arange(theta0, 2 * np.pi, delta_theta)

    thetas = thetas.reshape((len(thetas), 1))

    symbols = R_APSK * np.exp(1j * thetas)
    labels = grayCoded(ring_num) + ring_num * grayCoded(int(M_APSK[0])).T

    symbols = symbols.reshape((M, 1))
    labels = labels.reshape((M, 1))

    symbols = normalizeConste(symbols)

    if not dec_label:  # 转换为二进制标识
        labels = de2bi(labels, K)

    return symbols, labels


def grayConverter(label1, label2, K1, K2):

    label1_bin = de2bi(label1, K1)
    label2_bin = de2bi(label2, K2)

    label2_bin_converted = label2_bin

    for i in range(0, K1, 2):
        label2_bin_converted[0] = int(not (label2_bin_converted[0] ^ label1_bin[i]))

    for i in range(1, K1, 2):
        label2_bin_converted[1] = int(not (label2_bin_converted[1] ^ label1_bin[i]))

    return label2_bin_converted @ (2 ** np.array(range(K2 - 1, -1, -1)))


def MUST(symbols1, symbols2, labels1, labels2, P1, P2, MUST_type=1, dec_label=False):

    M1 = symbols1.size
    M2 = symbols2.size
    M = M1 * M2
    K1 = int(np.log2(M1))
    K2 = int(np.log2(M2))
    K = K1 + K2

    if MUST_type == 1:

        symbols = np.sqrt(P1) * symbols1 + np.sqrt(P2) * symbols2.T
        symbols = np.reshape(symbols, (M, 1))
        labels = labels1 * M2 + labels2.T
        labels = np.reshape(labels, (M, 1))

    elif MUST_type == 2:

        symbols = np.zeros((M, 1), dtype=complex)
        labels = np.zeros((M, 1), dtype=np.int32)

        I = labels2.argsort(axis=0)
        labels2 = labels2[I].reshape((M2, 1))
        symbols2 = symbols2[I].reshape((M2, 1))

        for i in range(M1):

            label1 = labels1[i]
            symbol1 = symbols1[i]

            for j in range(M2):

                label2 = labels2[j]
                labels[i * M2 + j] = label1 * M2 + label2
                label2_gray_converted = grayConverter(label1, label2, K1, K2)
                symbols[i * M2 + j] = (
                    symbol1 * np.sqrt(P1)
                    + np.sqrt(P2) * symbols2[label2_gray_converted]
                )

    else:
        pass

    if not dec_label:
        labels = de2bi(labels, K)

    return symbols, labels


def showConsteDiagram(symbols, labels, show_labels=True):

    M = symbols.size

    fig = plt.figure("Normalized Constellation Diagram", (10, 10))
    ax = axisartist.Subplot(fig, 1, 1, 1)
    fig.add_axes(ax)
    ax.axis[:].set_visible(False)
    ax.axis["x"] = ax.new_floating_axis(0, 0)
    ax.axis["y"] = ax.new_floating_axis(1, 0)
    ax.axis["x"].set_axis_direction("top")
    ax.axis["y"].set_axis_direction("left")
    ax.axis["x"].set_axisline_style("-|>", size=2.0)
    ax.axis["y"].set_axisline_style("-|>", size=2.0)
    ax.axis["x"].line.set_linewidth(2)
    ax.axis["y"].line.set_linewidth(2)
    ax.set_xticks([])
    ax.set_yticks([])

    label_font = {
        "family": "Times New Roman",
        "weight": "light",
        "size": 15,
    }

    for i in range(M):

        symbol, label = symbols[i], labels[i]

        label_string = ""  # 得到标识对应的字符串
        for i in range(len(label)):
            label_string += str(label[i])

        ax.scatter(symbol.real, symbol.imag, s=80, c="k")

        if show_labels:
            ax.text(
                symbol.real - 0.03,
                symbol.imag + 0.03,
                label_string,
                fontdict=label_font,
                color="k",
            )


def showSCDiagram(symbols, labels, M1, M2, show_labels=True):

    M = symbols.size

    fig = plt.figure("Normalized Constellation Diagram", (10, 10))
    ax = axisartist.Subplot(fig, 1, 1, 1)
    fig.add_axes(ax)
    ax.axis[:].set_visible(False)
    ax.axis["x"] = ax.new_floating_axis(0, 0)
    ax.axis["y"] = ax.new_floating_axis(1, 0)
    ax.axis["x"].set_axis_direction("top")
    ax.axis["y"].set_axis_direction("left")
    ax.axis["x"].set_axisline_style("-|>", size=2.0)
    ax.axis["y"].set_axisline_style("-|>", size=2.0)
    ax.axis["x"].line.set_linewidth(2)
    ax.axis["y"].line.set_linewidth(2)
    ax.set_xticks([])
    ax.set_yticks([])

    label_font = {
        "family": "Times New Roman",
        "weight": "light",
        "size": 15,
    }

    color_table = ["r", "g", "b", "c", "m", "y", "k"]
    col_idx = -1

    for i in range(M1):

        col_idx = (col_idx + 1) % len(color_table)
        # col = color_table[col_idx]
        col = "k"

        for j in range(M2):

            symbol, label = symbols[i * M2 + j], labels[i * M2 + j]

            label_string = ""  # 得到标识对应的字符串
            for k in range(len(label)):
                label_string += str(label[k])

            ax.scatter(symbol.real, symbol.imag, s=80, c=col)

            if show_labels:
                ax.text(
                    symbol.real - 0.3,
                    symbol.imag + 0.25,
                    label_string,
                    fontdict=label_font,
                    color=col,
                )


def getConsteQuad(symbols, labels, quad_num):

    if quad_num == 1:
        a, b = 1, 1
    elif quad_num == 2:
        a, b = -1, 1
    elif quad_num == 3:
        a, b = -1, -1
    else:
        a, b = 1, -1

    labels_quad = labels[a * symbols.real > 0]
    symbols_quad = symbols[a * symbols.real > 0]

    labels_quad = labels_quad[b * symbols_quad.imag > 0]
    symbols_quad = symbols_quad[b * symbols_quad.imag > 0]

    if symbols_quad.size > 1:
        symbols_quad -= symbols_quad.mean()

    symbols_quad = normalizeConste(symbols_quad)
    M2 = symbols_quad.size
    labels_quad = labels_quad % M2

    symbols_quad = symbols_quad.reshape((M2, 1))
    labels_quad = labels_quad.reshape((M2, 1))

    return symbols_quad, labels_quad


def getAdjustedConste(r, PHI):

    m = PHI.shape[0]
    M1 = 2 ** (m + 2)
    M2 = r.size

    M = M1 * M2

    _, labels = grayAPSK(np.array([M1] * M2))

    delta_theta = 2 * np.pi / M1
    theta0 = 0.5 * delta_theta

    thetas = np.arange(theta0, 2 * np.pi, delta_theta)
    thetas = thetas.reshape((len(thetas), 1))

    tmp1 = np.array([range(M1)]) + 1
    tmp2 = np.array([range(m)]) + 1

    phase_shifts = (-1) ** (np.ceil(tmp1.T / 2 ** (tmp2 - 1)) + 1) @ PHI

    THETAS = thetas + phase_shifts

    symbols = np.exp(1j * THETAS) * r
    symbols = np.reshape(symbols, (M, 1))

    symbols = normalizeConste(symbols)

    return symbols, labels


def decision(i, L, bit_capacity_1, bit_capacity_2, opt):

    metric = bit_capacity_2 / bit_capacity_1

    if opt == "good":  # good strategy
        assign_order = metric.argsort()
    elif opt == "bad":  # bad strategy
        assign_order = metric.argsort()[::-1]

    idx = i // L

    return assign_order[idx]


# 插值拟合计算信噪比门限
def getSNRth(SNR_dB, R, R_target):

    f = interpolate.interp1d(R, SNR_dB, kind="linear")

    SNRth = f(R_target)

    return SNRth.item()


def getBDMRate(bit_capacity, L, bit_num1, bit_num2, opt):
    # 在给定比特分配方案下，不同SNR下的两个用户的可达速率

    K = bit_capacity.shape[1]

    if opt == "good":
        assign_order = bit_capacity[0, :].argsort()[::-1]  # 从大到小
    elif opt == "bad":
        assign_order = bit_capacity[0, :].argsort()  # 从小到大
    else:
        print("未定义选项")
        return np.inf, np.inf

    if bit_num1 + bit_num2 != L * K:
        print("比特分配有问题，请检查")
        return np.inf, np.inf

    R_L = bit_capacity.sum(axis=1) * L

    d, m = divmod(bit_num2, L)
    R2_L = (
        bit_capacity[:, assign_order[range(d)]].sum(axis=1) * L
        + bit_capacity[:, assign_order[d]] * m
    )

    R1_L = R_L - R2_L

    R1 = R1_L / L
    R2 = R2_L / L

    return R1, R2


def getTheoreticRate(SNR1_dB, SNR2_dB, alpha):

    SNR1 = 10 ** (0.1 * SNR1_dB)
    SNR2 = 10 ** (0.1 * SNR2_dB)

    R1 = np.log2(1 + SNR1 * (1 - alpha) / (SNR1 * alpha + 1))
    R2 = np.log2(1 + SNR2 * alpha)

    return R1, R2


def getTheoreticSNRth(R1, R2, alpha):

    SNRth1 = (2 ** R1 - 1) / (1 - alpha)
    SNRth2 = (2 ** R2 - 1) / (1 - (1 - alpha) * (2 ** R2))

    return 10 * np.log10(SNRth1), 10 * np.log10(SNRth2)


def getMSDCMRatePair(symbols1, symbols2, SNR1_lin, SNR2_lin, alpha):

    M1 = symbols1.size
    M2 = symbols2.size

    M = M1 * M2

    device = torch.device("cuda:0")

    xx, ww = GH(n=16)
    xx = xx.reshape((xx.size, 1, 1))
    xx = torch.tensor(xx, dtype=torch.complex64, device=device)
    ww = torch.tensor(ww, dtype=torch.float32, device=device)

    # 基本层
    s1 = torch.tensor(
        symbols1 * np.sqrt(SNR1_lin * (1 - alpha)), dtype=torch.complex64, device=device
    )
    s2 = torch.tensor(
        symbols2 * np.sqrt(SNR1_lin * alpha), dtype=torch.complex64, device=device
    )

    symbols = s1 + s2.T
    symbols = symbols.reshape((M, 1))

    metric_num = torch.exp(-abs(xx - symbols + symbols.T) ** 2)
    metric_den = torch.exp(-abs(xx - s2 + s2.T) ** 2)

    tmp_num = torch.log2(torch.sum(metric_num, axis=1))
    tmp_den = torch.log2(torch.sum(metric_den, axis=1))

    tmp1 = torch.sum(ww * torch.sum(tmp_num, axis=1)) / (M * np.pi)
    tmp2 = torch.sum(ww * torch.sum(tmp_den, axis=1)) / (M2 * np.pi)

    Rbc = np.log2(M1) - tmp1 + tmp2

    # 增强层
    s2 = torch.tensor(
        symbols2 * np.sqrt(SNR2_lin * alpha), dtype=torch.complex64, device=device
    )

    metric2 = torch.exp(-abs(xx - s2 + s2.T) ** 2 + abs(xx) ** 2)

    tmp = torch.log2(torch.sum(metric2, axis=1))

    Ruc = np.log2(M2) - torch.sum(ww * torch.sum(tmp, axis=1)) / (M2 * np.pi)

    return Rbc, Ruc


def getMSDBICMRatePair(symbols1, labels1, symbols2, labels2, SNR1_lin, SNR2_lin, alpha):

    M1 = symbols1.size
    M2 = symbols2.size
    K1 = int(np.log2(M1))
    K2 = int(np.log2(M2))

    K = K1 + K2
    M = M1 * M2

    device = torch.device("cuda:0")

    xx, ww = GH(n=16)
    xx = xx.reshape((xx.size, 1, 1))
    xx = torch.tensor(xx, dtype=torch.complex64, device=device)
    ww = torch.tensor(ww, dtype=torch.float32, device=device)

    # 基本层容量推导
    symbols = np.sqrt(SNR1_lin) * (
        symbols1 * np.sqrt(1 - alpha) + symbols2.T * np.sqrt(alpha)
    )
    symbols = np.reshape(symbols, (M, 1), order="F")

    symbols = torch.tensor(symbols, dtype=torch.complex64, device=device)

    labels = labels1 * M2 + labels2.T
    labels = np.reshape(labels, (M, 1), order="F")
    labels = de2bi(labels, K)

    metric = torch.exp(-abs(xx - symbols + symbols.T) ** 2)

    Rbc = 0

    for k in range(K1):

        mask0 = np.array(np.ix_(labels[:, k] == 0))
        mask1 = np.array(np.ix_(labels[:, k] == 1))

        p00 = torch.sum(metric[:, mask0.T, mask0], 1)
        p10 = torch.sum(metric[:, mask1.T, mask0], 1)
        p01 = torch.sum(metric[:, mask0.T, mask1], 1)
        p11 = torch.sum(metric[:, mask1.T, mask1], 1)

        Rbc += 1 - torch.sum(
            ww
            * torch.sum(torch.log2(1 + p10 / p00) + torch.log2(1 + p01 / p11), axis=1)
        ) / (M * np.pi)

    del metric

    # 增强层
    s2 = np.sqrt(SNR2_lin) * torch.tensor(
        symbols2 * np.sqrt(alpha), dtype=torch.complex64, device=device
    )
    labels2 = de2bi(labels2, K2)

    metric = torch.exp(-abs(xx - s2 + s2.T) ** 2)

    Ruc = 0

    for k in range(K2):

        mask0 = np.array(np.ix_(labels2[:, k] == 0))
        mask1 = np.array(np.ix_(labels2[:, k] == 1))

        p00 = torch.sum(metric[:, mask0.T, mask0], 1)
        p10 = torch.sum(metric[:, mask1.T, mask0], 1)
        p01 = torch.sum(metric[:, mask0.T, mask1], 1)
        p11 = torch.sum(metric[:, mask1.T, mask1], 1)

        Ruc += 1 - torch.sum(
            ww
            * torch.sum(torch.log2(1 + p10 / p00) + torch.log2(1 + p01 / p11), axis=1)
        ) / (M2 * np.pi)

    return Rbc, Ruc


def awgn(symbols, labels, Kb=0):

    M = symbols.size
    K = int(np.log2(M))

    if Kb == 0:
        Kb = K
        Krange = range(Kb)
    elif Kb > 0:
        Krange = range(Kb)
    else:
        Kb = -Kb
        Krange = range(K - Kb, K)

    if labels.shape[1] != K:
        labels = de2bi(labels, K)

    device = torch.device("cuda:0")

    xx, ww = GH(n=16)
    xx = xx.reshape((xx.size, 1, 1))
    xx = torch.tensor(xx, dtype=torch.complex64, device=device)
    ww = torch.tensor(ww, dtype=torch.float32, device=device)

    symbols = torch.tensor(symbols, dtype=torch.complex64, device=device)

    metric = torch.exp(-abs(xx - symbols + symbols.T) ** 2)

    bit_capacity = np.zeros(Kb)

    for i in range(len(Krange)):

        k = Krange[i]
        mask0 = np.array(np.ix_(labels[:, k] == 0))
        mask1 = np.array(np.ix_(labels[:, k] == 1))

        p00 = torch.sum(metric[:, mask0.T, mask0], 1)
        p10 = torch.sum(metric[:, mask1.T, mask0], 1)
        p01 = torch.sum(metric[:, mask0.T, mask1], 1)
        p11 = torch.sum(metric[:, mask1.T, mask1], 1)

        bit_capacity[i] = 1 - torch.sum(
            ww
            * torch.sum(torch.log2(1 + p10 / p00) + torch.log2(1 + p01 / p11), axis=1)
        ) / (M * np.pi)

    return bit_capacity


def de2bi(nums, K):

    if nums.size > 1:
        nums = nums.squeeze()  # 将nums尺寸统一为(M, )

    M = len(nums)

    ret = []
    for i in range(M):  # 遍历所有要转换的十进制数
        num = nums[i]
        tmp = []
        for i in range(K - 1, -1, -1):  # 逐二进制位遍历
            tmp.append((num >> i) & 1)
        ret.append(tmp)

    ret = np.array(ret)

    if ret.shape[0] == 1:
        ret = ret.squeeze()

    return ret


def GH(n=16):

    p = np.zeros((n + 1, n + 1))
    p[0, n] = 1
    p[1, n - 1 :] = np.array([1, 0]) * 2
    for k in range(2, n + 1):
        p[k, n - k :] = 2 * np.append(p[k - 1, n - k + 1 :], 0) - 2 * (
            k - 1
        ) * np.insert(p[k - 2, n - k + 2 :], 0, [0, 0])

    x = np.roots(p[n, :])
    w = np.zeros(n)
    for k in range(n):
        w[k] = (
            2 ** (n - 1)
            * np.math.factorial(n)
            * np.sqrt(np.pi)
            / (n ** 2 * np.polyval(p[n - 1, :], x[k]) ** 2)
        )

    x = np.reshape(x, (1, n))
    w = np.reshape(w, (1, n))

    xx = x + 1j * x.T
    ww = w * w.T

    xx = xx.flatten(order="F")
    ww = ww.flatten(order="F")

    xx = xx.reshape((1, len(xx)))
    ww = ww.reshape((1, len(ww)))

    return xx, ww
