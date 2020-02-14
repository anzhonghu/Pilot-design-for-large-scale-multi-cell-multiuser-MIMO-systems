M=[1;5;11;17;23;25;31];
pilots_all2 = zeros(N_p,K*L);
pilots = zeros(N_p,K);%Chu sequences
for i = 1 : L
    for k = 1 : K
        for l = 0 : N_p-1
            ll = mod(l+k-1, N_p);
            pilots(l+1,k) = exp(1i*pi*M(i,1)*ll*ll/N_p);
        end
    end
    pilots_all2(:,(i-1)*K+1:i*K) = pilots;
end