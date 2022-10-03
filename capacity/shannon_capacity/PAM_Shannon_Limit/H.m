function entropy = H(p)

    if p<1e-10
        p=1e-10;
    end
    if length(p)==1
        entropy = -p*log2(p)-(1-p)*log2(1-p);
    else
        entropy = -sum(p.*log2(p));
    end
end

