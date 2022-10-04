function capacity = SISO_Capacity(conste, SNR, chan_type, CM_type)

    symbols = conste.symbols;
    M = length(symbols);

    if isfield(conste, 'prob')
        syms_prob = conste.prob;
    else
        syms_prob = ones(M, 1) / M;
    end

    if CM_type == 'BICM'

        if isfield(conste, 'labels')
            labels = conste.labels;
        else
            error('Labels of constellation points must be specified!');
        end

    end

    if any(imag(symbols))
        conste_dim = 2;
    else
        conste_dim = 1;
    end

    SNR_type = SNR.SNR_type;
    SNR_vec = SNR.SNR_vec;
