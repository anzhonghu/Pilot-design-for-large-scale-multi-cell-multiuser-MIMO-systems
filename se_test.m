clear;
close all;
%%constants
Ns = 14;
L = 7;%cell number
N_p = 3 * Ns;%pilot length
K = N_p - 15;%user number,27
ge_sequence1;%generate designed pilots
bound = (K*L)^2/N_p - K*L;%generate welch bound
bcor = 0;
for i = 1 : L
    for j = 1 : L
        if j ~= i
            A = pilots_all(:,(j-1)*K+1:j*K).' * conj(pilots_all(:,(i-1)*K+1:i*K)) / N_p;
            X = A.*conj(A);
            bcor = bcor + sum(sum(X));
        else
        end
    end
end
