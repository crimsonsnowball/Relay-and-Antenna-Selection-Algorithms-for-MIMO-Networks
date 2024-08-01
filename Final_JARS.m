clear;clc;

%% Analysis

for K=1:4
    i = 1;
    for E = 5:1:25
        F_T1(i) = (1-exp(-15/(10.^(E/10))).*(1 + 15/(10.^(E/10)) ))^2;
        F_T2(i) = (1-exp(-15/(10.^(E/10)))*(1+15/(10.^(E/10))))^2;
        i=i+1; % R_bar = 15 , Ms = MRk = 2
    end
    y(K,:)= (F_T1 + F_T2 - F_T2.*F_T1).^K;
end

E = 5:1:25;
semilogy(E,y)
ylim([1e-6 1])
grid on
hold on

%% Simulation

MRk = 2; % Total Relay antenna
MS = 2; % Total antenna at Source
MD = 2; % Total antenna at Destination
R = 2; % e2e Spectral Efficency

N = 10^6; % Ittrations

for K = 1:4

    q=1;
    for P_source = 5:1:25
    
        % Source to Relay
        h_sr_real = normrnd(0,1/sqrt(2),MS*K,MRk,N);
        h_sr_img = normrnd(0,1/sqrt(2),MS*K,MRk,N);
        C1 = h_sr_real + 1i*h_sr_img; % C1 = Channel Matrix form Source to Relay 
        gamma_ms_mrk = 10^(P_source/10).*(abs(C1)).^2;  % P*|h|^2 where P is in Watts 
        gamma_ms = sum(gamma_ms_mrk,2); % MRC at Relay = Sum of all SNRs
        
        % Realy to Destination
        h_rd_real = normrnd(0,1/sqrt(2),MRk*K,MD,N);
        h_rd_img = normrnd(0,1/sqrt(2),MRk*K,MD,N);
        C2 = h_rd_real + 1i*h_rd_img; % h ~ CN(0,1) Source to Relay Channel 
        gamma_mrk_md = 10^(P_source/10)*(abs(C2)).^2;  % P*|h|^2 where P is in Watts
        gamma_mrk = sum(gamma_mrk_md,2); % MRC at Destination = Sum of SNRs
        
        j = 1;
        
        for i=1:MS:2*K % one loop for Kth relay
            gamma_ms_Rk(:,:,:) = [gamma_ms(i,:,:);gamma_ms(i+1,:,:)];
            gamma_mrk_Rk(:,:,:) = [gamma_mrk(i,:,:);gamma_mrk(i+1,:,:)]; % selecting two gamma at a time
            
            Z_Rk = [repmat(gamma_ms_Rk,2,1),repelem(gamma_mrk_Rk,2,1)]; % All possible permutations of gamma(ms) and gamma(mrk)
            
            gamma_e2e_min_Rk = min(Z_Rk,[],2); % min of Z row wise
            
            gamma_e2e_max_Rk_N(j,:) = max(gamma_e2e_min_Rk); % max % Column = N=1,2,3.. , Row = R1,R2,R3,R4
            j=j+1;
        end
        
        gamma_max_outofall_relays = max(gamma_e2e_max_Rk_N,[],1); % Best SNR in Nth ittration
        
        Outage_logic_array = gamma_max_outofall_relays < 2^(2*R) -1; % checks if outage has occured in Nth itteration
        Prob_outage(K,q) = sum(Outage_logic_array)/N; % Outage prob of system for Qth value of Source Power
        q=q+1;
    end

end

P_source = 5:1:25;
semilogy(P_source,Prob_outage,'diamond')
ylim([1e-6 1])
xlabel('Source Transmit Power in dB');
ylabel('Outage Probability');
legend('K=1','K=2','K=3','K=4');
grid on
hold off