%Pilot design for multicell frequence reuse
%This file is for the simulation calculation
clear;
close all;
%%constants
Ns = 14;
L = 7;%cell number
N_p = 3 * Ns;%pilot length
K = N_p - 15;%user number,27
rc = 1600;%radius of the cell(center to edge)(m)
rh = 1500;%minimum terminal radius of the cell(m)
rcx = rc * 0.01;
rhx = rh * 0.01;
gamma = 3.8;%decay exponent
sigma = 8;%in dB
Num = 5e3;%Iteration number
% Num = 1e2;%Iteration number
ant_num = 5e2;
ge_sequence;%generate designed pilots
ge_sequence2;%generate  pilots in Kangguixia

%%position of every base
base(1:7,1) = [0;(1i * 2 * rc);(sqrt(3) * rc + 1i * rc);(sqrt(3) * rc - 1i * rc);(-1i * 2 * rc);(-sqrt(3) * rc - 1i * rc);(-sqrt(3) * rc + 1i * rc);];
SIR_pick = zeros(10,4);
count = 0;
for i_ant = 50 : 50 : ant_num
    count = count + 1;
    G = zeros(i_ant,K*L*L);
    for jj = 1 : Num
        ge_sequence1;
        dis(1:K,1:7) = (rem(rand(K,7),rcx-rhx) + rhx) * 100;%, if rh=100, the intracell interference is stronger than inter-cell interference,
        % so the pilots of kang's will be better.  Hence, rh=1500 here.
        ang(1:K,1:7) = rand(K,7) * 2 * pi;
        pos(1:K,1:7) = dis .* (exp(1i * ang));
        pos(:,2) = pos(:,2) + base(2,1);
        pos(:,3) = pos(:,3) + base(3,1);
        pos(:,4) = pos(:,4) + base(4,1);
        pos(:,5) = pos(:,5) + base(5,1);
        pos(:,6) = pos(:,6) + base(6,1);
        pos(:,7) = pos(:,7) + base(7,1);
        SIR = zeros(K*7,4);
        shadow_amp1 = sqrt(10.^(randn(1,K*L) * sigma * 0.1));
%         shadow_amp1 = ones(i_ant,K*L);
        %%Interference for 1st /3rd kind
        for l1 = 1 : L%BS
            for l2 = 1 : L%user
                G(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = 1 / sqrt(2) * (randn(i_ant,K)+1i*randn(i_ant,K))...
                    * sqrt(diag(shadow_amp1(:,(l2-1)*K+1:l2*K))  / diag((abs(pos(:,l2)-base(l1,1))*0.01).^gamma));
            end
        end
        a = sign(randn(K,L));%data
        for k = 1 : K
            for l = 1 : L
                if 0 == a(k,l)
                    a(k,l) = 1;
                else
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %as the attenuation between cells is large, the channel matrices are
        %not quite orthogonal, some terminals are quite close to the cell,
        %so the intercell interference is smaller than intracell interference.
        for j = 1 : L%BS
            Gjj = G(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K);
            xla =  zeros(i_ant,1);%data
            xl1 =  zeros(i_ant,K);
            xl2 =  zeros(i_ant,K);
            xl3 =  zeros(i_ant,K);
            xl4 =  zeros(i_ant,K);
            for l = 1 : L
                ll = l;
                jn = ll;
                Gjl = G(:,(j-1)*L*K+(l-1)*K+1:(j-1)*L*K+l*K);
                xla = xla + Gjl * a(:,l);%data
                xl1 = xl1 + Gjl * pilots_all(:,(j-1)*K+1:j*K).' * conj(pilots_all(:,(l-1)*K+1:l*K)) / N_p;
                xl2 = xl2 + Gjl;
                xl3 = xl3 + Gjl * pilots_all1(:,(j-1)*K+1:j*K).' * conj(pilots_all1(:,(l-1)*K+1:l*K)) / N_p;
                xl4 = xl4 + Gjl * pilots_all2(:,(j-1)*K+1:j*K).' * conj(pilots_all2(:,(l-1)*K+1:l*K)) / N_p;
            end
            cor_Gjj = 1/i_ant * (Gjj' * Gjj);
            z(:,1) = 1/i_ant * (xl1)' * (xla);%proposed
            z(:,2) = 1/i_ant * (xl2)' * (xla);%1st, original,and is equivalent to the chu proposed for vehicular environments by Wangwenjun.
            z(:,3) = 1/i_ant * (xl3)' * (xla);%random
            z(:,4) = 1/i_ant * (xl4)' * (xla);%Kangguixia
            x = zeros(K,1);
            for k = 1 : K
                x(k,1) =  cor_Gjj(k,k) * a(k,j);%%This cycle is particularly important,otherwise, the SIR will be wrong
            end
            n(:,1) = z(:,1) - x;%the inter-cell interference is too small, so the pilots have little effects
            n(:,2) = z(:,2) - x;
            n(:,3) = z(:,3) - x;
            n(:,4) = z(:,4) - x;
            for k = 1 : K
                SIR((j-1)*K+k,1) = (abs(x(k,1)))^2 / (abs(n(k,1)))^2;
                SIR((j-1)*K+k,2) = (abs(x(k,1)))^2 / (abs(n(k,2)))^2;
                SIR((j-1)*K+k,3) = (abs(x(k,1)))^2 / (abs(n(k,3)))^2;
                SIR((j-1)*K+k,4) = (abs(x(k,1)))^2 / (abs(n(k,4)))^2;
            end
        end
        %%combine all
        SIR_s((jj-1)*K*7+1:jj*K*7,:) = SIR;
    end
    SIR_s = 10 * log10(SIR_s);
    SIR_s = sort(SIR_s);
    SIR_pick(count,:) = SIR_s(Num*K*7*0.05,:);
    sprintf('%d',i_ant)
end

ant_num=5e2;
antenna = (50:50:ant_num);
h = figure;
xx = axes('FontSize',16);
set(h,'PaperType','A4');
plot(antenna,SIR_pick(:,4),'k-.x','LineWidth',2,'MarkerSize',10)
hold on
plot(antenna,SIR_pick(:,2),'k-.s','LineWidth',2,'MarkerSize',10)
plot(antenna,SIR_pick(:,3),'k-.o','LineWidth',2,'MarkerSize',10)
plot(antenna,SIR_pick(:,1),'k-+','LineWidth',2,'MarkerSize',10)
le = legend('Pilots in [8]','Pilots in [9]','Random pilots','Proposed pilots','Location','SouthEast');
set(le,'Fontsize',16,'Fontname','Times')
xlabel('Number of BS antennas','Fontsize',20,'Fontname','Times')
ylabel('P \geq 0.95 SIR (dB)','Fontsize',20,'Fontname','Times')
%print(h,'-dpdf','7_cell_SIR_antenna')

