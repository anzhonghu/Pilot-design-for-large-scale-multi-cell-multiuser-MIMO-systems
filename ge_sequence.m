
c = N_p - K;
delta = c - L + 4;%corresponds to c_d in paper, so +2 is changed to +4

pilots = zeros(N_p,K);%Chu sequences
for k = 1 : K
    for l = 0 : N_p-1
        ll = mod(l+k-1, N_p);
        pilots(l+1,k) = exp(1i*pi*ll*ll/N_p);
    end
end
R_p = pilots.' * conj(pilots) / N_p;

b1 = zeros(N_p,K);
for k = 1 : K
    for l = 0 : N_p-1
        ll = mod(l+k-1, N_p);
        b1(l+1,k) = pilots(l+1,k) * exp(1i*2*pi*ll/N_p);
    end
end
R_b1p = b1.' * conj(pilots) / N_p;

b2 = zeros(N_p,K);
for k = 1 : K
    for l = 0 : N_p-1
        ll = mod(l+k-1, N_p);
        b2(l+1,k) = pilots(l+1,k) * exp(1i*2*pi*2*ll/N_p);
    end
end
R_b2p = b2.' * conj(pilots) / N_p;

b3 = zeros(N_p,K);
for k = 1 : K
    for l = 0 : N_p-1
        ll = mod(l+k-1, N_p);
        b3(l+1,k) = pilots(l+1,k) * exp(1i*2*pi*delta*ll/N_p);
    end
end
R_b3p = b3.' * conj(pilots) / N_p;

b4 = zeros(N_p,K);
for k = 1 : K
    for l = 0 : N_p-1
        ll = mod(l+k-1, N_p);
        b4(l+1,k) = pilots(l+1,k) * exp(1i*2*pi*(delta+1)*ll/N_p);
    end
end
R_b4p = b4.' * conj(pilots) / N_p;

b5 = zeros(N_p,K);
for k = 1 : K
    for l = 0 : N_p-1
        ll = mod(l+k-1, N_p);
        b5(l+1,k) = pilots(l+1,k) * exp(1i*2*pi*(delta+2)*ll/N_p);
    end
end
R_b5p = b5.' * conj(pilots) / N_p;

b6 = zeros(N_p,K);
for k = 1 : K
    for l = 0 : N_p-1
        ll = mod(l+k-1, N_p);
        b6(l+1,k) = pilots(l+1,k) * exp(1i*2*pi*(delta+3)*ll/N_p);
    end
end
R_b6p = b6.' * conj(pilots) / N_p;

pilots_all = [pilots,b1,b2,b3,b4,b5,b6];
