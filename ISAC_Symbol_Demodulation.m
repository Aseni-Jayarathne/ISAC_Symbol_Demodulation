clear
clc;
close all;

% Parameters
L = 1e5; %data size
M = 4; %modulation order

NTx = 4; %number of transmit antennas
NRx = [1,2,4]; %number of receiver antennas

phase_shift = 0; %artificial phase shift
constellation = pskmod(0:(M-1), M, phase_shift); %Alphabet
beta = logspace(-3, 3, 20); %power of communication signal
gamma = 1; %power of sensing signal
theta = 30; %target angle
theta_rad = deg2rad(theta); %target angle in rad
d = 1/2; % Spacing between antennas (in wavelength units)
sigma_sq_dB = -30; %noise varience in dB
sigma_sq = 10.^(sigma_sq_dB./10); %noise varience
p = 1; %number of minimum points for Maximum Likelihood detection



%Data

%UE
data_UE = randi([0 M-1],1, L);
data_bin_UE = de2bi(data_UE,log2(M)); 

txSig_UE = pskmod(data_UE,M,phase_shift); %Transmitted signal from UE


% Target
data_T = randi([0 1],1, L);
data_bin_T = de2bi(data_T,log2(M));

txSig_T = pskmod(data_T,M,phase_shift); %Transmitted signal from target



for nrx = 1:length(NRx)
   
    
    % Beamforming
    alpha = sqrt(0.5)*(randn(1, L)+1i*randn(1, L)); %reflection coefficient vector
    a = exp(1i * 2 * pi * (0:NTx-1)* sin(theta_rad)*d).'/sqrt(NTx);% transmit steering vector
    b = exp(1i * 2 * pi * (0:NRx(nrx)-1)* sin(theta_rad)*d).'/sqrt(NRx(nrx)); % receive steering vector
    f = a; %transmit beamforming vector

    
    
    % Covarience between alpha and noise
    g = b*a'*f; %Sensing channel
    R = g*g' + sigma_sq*eye(NRx(nrx)); % covarience matrix between alpha and noise

    gg = g'*inv(R);

    % Independent subchannels
    Q = sqrtm(inv(R)); % Whitening matrix

    for n = 1:length(N)

        % Matrices to store data
        discrepancy = [];
        ber = [];
        mse_IC = [];
        mse_NIC = [];

 

        for ii=1:length(beta) %iterate through beta values
            fprintf('%d\n',ii);
     
    
            %User Equipment to BS
            rng(1)           
            
            h_c = sqrt(beta(ii))*sqrt(1/2).*(randn(NRx(nrx), L)+1i*randn(NRx(nrx),L)); %Rayleigh fading channel for LOS 


            % Maximum Likelihood Detection
            ML_predict_txSig = [];
            
            % Estimated alpha
            alpha_hat = zeros(1, L);
            alpha_mmse = zeros(1,L);
            tot_neu = 0;
            tot_den = 0;
            for jj = 1:L
                %Additive White Gaussian Noise
                noise = sqrt(sigma_sq/2) * (randn(NRx(nrx), 1) + 1i*randn(NRx(nrx), 1)); %AWGN noise with zero mean and sigma_sq variance
          
                %Received Data
                %Communication signal
                Comm_Rx = h_c(:,jj).*txSig_UE(:,jj); % received communication data
    
                %Sensing signal
                Sensing_Rx = sqrt(gamma)*alpha(:,jj).*(b*a'*f*txSig_T(:,jj)); % received sensing data
               
    
                rxSig = Comm_Rx + Sensing_Rx + noise; % Received signal vector befor combining
                 
    
                w = Q*h_c(:,jj); %MRC
                rxSig_MRC = w'*Q*rxSig;
                euclidean_dist = []; % Euclidean distance vector
    
                for k = 1:M %loop through constellation
                    ref= w'*Q*h_c(:,jj)*constellation(k);  
                    euclidean_dist = [euclidean_dist,norm(rxSig_MRC-ref)];
                end      
            
                [x_val, x_indx] = sort(euclidean_dist); %Minimum Euclidean distance
                predicted_symbol = x_indx(1:p) - 1;
                ML_predict_txSig(jj) = predicted_symbol; %Prediction
            

                % Sensing Performance
                s_hat = pskmod(predicted_symbol,M,phase_shift); %Transmitted signal from UE
                
                % Approximated Communication signal
                Comm_approx = h_c(:,jj).*s_hat; % received communication data
                y_hat_s = rxSig - Comm_approx;
                              
    
                ka = b*a'*f*txSig_T(:,jj); %Calculate ka matrix
                K = ka'*inv(ka'*ka + sigma_sq*eye(1));
                alpha_hat(jj) = K*y_hat_s;

                
            end


            % Compute Simulated BER
            predict_bin = de2bi(ML_predict_txSig,log2(M));
            discrepancy = [discrepancy,mean2(predict_bin ~= data_bin_UE)]; 
    
            %MSE of alpha
            mse_IC = [mse_IC, immse( alpha , alpha_hat )];

        end 
            
    
        
        %Plot Curves
        figure(1)
        loglog(beta,discrepancy,'-s','LineWidth',1.5,'DisplayName',strcat('N : ',num2str(N(n))));
        hold on
        grid on
        box on
        xlim([beta(1),beta(length(beta))]);
        legend('show','FontSize',14)
        xlabel("\beta (Communication Signal Power)")
        ylabel("BER",'FontSize')
        title ('SNR vs BER curve for a ISAC+RIS System')

    
        figure(2)
        loglog(beta,mse_IC,'-s','LineWidth',1.5, 'DisplayName',strcat('N: ',num2str(N(n))));
        hold on
        grid on
        box on
        legend('show','FontSize',14)
        xlim([beta(1),beta(length(beta))]);

        xlabel("\beta (Communication Signal Power)")
        ylabel("MSE of \alpha",'FontSize')
        title ('Sensing performance curve for a ISAC+RIS System')

    end
    


       
           

    


    
end

set(gca,'FontSize',15)





