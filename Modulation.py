import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist


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


def normalizeConste(symbols):

    M = symbols.size

    normalized_symbols = symbols / np.linalg.norm(symbols) * np.sqrt(M)

    return normalized_symbols


def PSK(M, labeling="gray"):

    symbols = np.exp(1j * 2 * np.pi / M * np.arange(M)).reshape(1, M)
    labels = grayCoded(M).reshape(1, M)

    if labeling == "gray":
        pass
    elif labeling == "bin":
        labels = np.arange(M).reshape(1, M)
    else:
        print("Unsupported labeling name! Gray labeling is used instead.")
        print("But you can customize the labeling with a vector outside this function.")

    return symbols, labels


def QAM(M, labeling="gray"):

    MI, MQ = M[0], M[1]

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

    return symbols, labels


def PAM(M, labeling="gray"):

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

    return symbols, labels


def showConsteDiagram(symbols, labels, show_labels=True):

    M = symbols.size

    labels = de2bi(labels)

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

