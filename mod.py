import numpy as np
import matplotlib.pyplot as plt 
import mpl_toolkits.axisartist as axisartist


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


def grayCoded(n):

    ret = np.zeros((1, n), dtype=np.int32)

    for i in range(n):
        ret[0, i] = i ^ (i >> 1)

    return ret


def normalizeConste(symbols):
    
    M = len(symbols)
    normalized_symbols = symbols/np.linalg.norm(symbols)*np.sqrt(M)

    return normalized_symbols


def adjustBinLabel(symbols, labels):

    M = symbols.size
    K = int(np.log2(M))
    
    tmp1 = symbols.real>0
    tmp2 = symbols.imag>0
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
    I_symbols = np.array([np.arange(-MI+1, MI+1, 2)])
    Q_symbols = -1j*I_symbols.T

    symbols = I_symbols+Q_symbols

    I_labels = grayCoded(MI) #维度为(1, MI)
    Q_labels = I_labels.T*MI #维度为(MI, 1)

    labels = I_labels+Q_labels #维度为(MI, MI)

    symbols = symbols.reshape((M, 1))
    labels = labels.reshape((M, 1))

    symbols = normalizeConste(symbols) #归一化能量
    
    labels = de2bi(labels, K)
    
    labels = adjustBinLabel(symbols, labels) #调整为两个高位比特决定象限
    
    if dec_label: #转换为二进制标识
        labels = labels@(2**np.array(range(K-1, -1, -1)))
        labels = labels.reshape((M, 1))
    
    return symbols, labels


def getAPSKRadius(m2):

    R_APSK = np.zeros((1, int(2**m2)))
    for i in range(R_APSK.size):
        R_APSK[0, i] = np.sqrt(-np.log(1-(i+1-0.5)*2**(-m2)))
        
    return R_APSK


def grayAPSK(M_APSK, dec_label=False):

    M_APSK = M_APSK.squeeze()
    
    M = np.sum(M_APSK) # 调制阶数
    K = int(np.log2(M))
    
    ring_num = M_APSK.size # APSK环数
    m2 = int(np.log2(ring_num))

    # 参见参考文献中 Fig.1
    theta0 = np.pi/M_APSK[0]
    delta_theta = 2*theta0

    R_APSK = getAPSKRadius(m2)
    thetas = np.arange(theta0, 2*np.pi, delta_theta)

    thetas = thetas.reshape((len(thetas), 1))

    symbols = R_APSK*np.exp(1j*thetas)
    labels = grayCoded(ring_num)+ring_num*grayCoded(int(M_APSK[0])).T

    symbols = symbols.reshape((M, 1))
    labels = labels.reshape((M, 1))

    symbols = normalizeConste(symbols)
    
    if not dec_label: #转换为二进制标识
        labels = de2bi(labels, K)
    
    return symbols, labels

def grayConverter(label1, label2, K1, K2):
    
    label1_bin = de2bi(label1, K1)
    label2_bin = de2bi(label2, K2)
   
    label2_bin_converted = label2_bin
    
    for i in range(0, K1, 2):
        label2_bin_converted[0] = int(not (label2_bin_converted[0]^label1_bin[i]))
   
    for i in range(1, K1, 2):
        label2_bin_converted[1] = int(not (label2_bin_converted[1]^label1_bin[i]))
       
    return label2_bin_converted@(2**np.array(range(K2-1, -1, -1)))
       

def MUST(symbols1, symbols2, labels1, labels2, P1, P2, MUST_type=1, dec_label=False):
    
    M1 = symbols1.size
    M2 = symbols2.size
    M = M1*M2
    K1 = int(np.log2(M1))
    K2 = int(np.log2(M2))
    K = K1+K2
    
    if MUST_type == 1:
        
        symbols = np.sqrt(P1)*symbols1+np.sqrt(P2)*symbols2.T
        symbols = np.reshape(symbols, (M, 1))
        labels = labels1*M2+labels2.T
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
                labels[i*M2+j] = label1*M2+label2
                label2_gray_converted = grayConverter(label1, label2, K1, K2)
                symbols[i*M2+j] = symbol1*np.sqrt(P1)+np.sqrt(P2)*symbols2[label2_gray_converted]
                
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
    ax.axis["x"].set_axis_direction('top') 
    ax.axis["y"].set_axis_direction('left') 
    ax.axis["x"].set_axisline_style("-|>", size = 2.0) 
    ax.axis["y"].set_axisline_style("-|>", size = 2.0) 
    ax.axis['x'].line.set_linewidth(2)
    ax.axis['y'].line.set_linewidth(2)
    ax.set_xticks([]) 
    ax.set_yticks([]) 
    
    label_font = {'family': 'Times New Roman',
             'weight': 'light',
             'size': 15,
             }
    
    for i in range(M):
        
        symbol, label = symbols[i], labels[i]
        
        label_string = '' # 得到标识对应的字符串
        for i in range(len(label)):
            label_string += str(label[i])
            
        ax.scatter(symbol.real, symbol.imag, s=80, c='k')
        
        if show_labels:
            ax.text(symbol.real-0.03, symbol.imag+0.03, 
                    label_string, fontdict=label_font, color='k')

def showSCDiagram(symbols, labels, M1, M2, show_labels=True):

    M = symbols.size

    fig = plt.figure("Normalized Constellation Diagram", (10, 10)) 
    ax = axisartist.Subplot(fig, 1, 1, 1) 
    fig.add_axes(ax) 
    ax.axis[:].set_visible(False) 
    ax.axis["x"] = ax.new_floating_axis(0, 0) 
    ax.axis["y"] = ax.new_floating_axis(1, 0) 
    ax.axis["x"].set_axis_direction('top') 
    ax.axis["y"].set_axis_direction('left') 
    ax.axis["x"].set_axisline_style("-|>", size = 2.0) 
    ax.axis["y"].set_axisline_style("-|>", size = 2.0) 
    ax.axis['x'].line.set_linewidth(2)
    ax.axis['y'].line.set_linewidth(2)
    ax.set_xticks([]) 
    ax.set_yticks([]) 
    
    label_font = {'family': 'Times New Roman',
             'weight': 'light',
             'size': 15,
             }
    
    color_table = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    col_idx = -1
    
    for i in range(M1):
        
        col_idx = (col_idx+1)%len(color_table)
        # col = color_table[col_idx]
        col = 'k'
        
        for j in range(M2):
        
            symbol, label = symbols[i*M2+j], labels[i*M2+j]
            
            label_string = '' # 得到标识对应的字符串
            for k in range(len(label)):
                label_string += str(label[k])
                
            ax.scatter(symbol.real, symbol.imag, s=80, c=col)
            
            if show_labels:
                ax.text(symbol.real-0.3, symbol.imag+0.25, 
                        label_string, fontdict=label_font, color=col)

def getConsteQuad(symbols, labels, quad_num):

    if quad_num==1:
        a, b = 1, 1
    elif quad_num==2:
        a, b = -1, 1
    elif quad_num==3:
        a, b = -1, -1
    else:
        a, b = 1, -1

    labels_quad = labels[a*symbols.real>0]
    symbols_quad = symbols[a*symbols.real>0]

    labels_quad = labels_quad[b*symbols_quad.imag>0]
    symbols_quad = symbols_quad[b*symbols_quad.imag>0]

    if symbols_quad.size>1:
        symbols_quad -= symbols_quad.mean()
        
    symbols_quad = normalizeConste(symbols_quad)    
    M2 = symbols_quad.size
    labels_quad = labels_quad%M2
    
    symbols_quad = symbols_quad.reshape((M2, 1))
    labels_quad = labels_quad.reshape((M2, 1))

    return symbols_quad, labels_quad