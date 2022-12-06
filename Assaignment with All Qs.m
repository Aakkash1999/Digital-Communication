clear all;
clc ;

%% generating a random BPSK stream of length 1000
L = 1000; % length of the bit sequnece
fs = 10; % sampling frequency
Tb = 1/ fs ;
dt = 0.001;
t = 0: dt : L / fs ; % time vector
bit_seq = randi ([0 1] ,1 , L ) ; % original bit sequence

%% mapping the bit stream into BPSK symbols
bpsk_sym = zeros (1 , L ) ;

for i = 1: L
    if ( bit_seq ( i ) ==0)
        bpsk_sym ( i ) = -1;
    else
        bpsk_sym ( i ) = 1;
    end
end

% plotting the bit sequence and the BPSK symbol train
bit_mod = zeros (1 , length ( t ) ) ;
for j = 0: L -1
    bit_mod ( j * Tb / dt +1:( j +1) * Tb / dt ) = bit_seq ( j +1) ;
end

BPSK_time = zeros (1 , length ( t ) ) ;
for n = 0: L -1
    BPSK_time ( round ( n * Tb / dt ) +1) = bpsk_sym ( n +1) ;
end

figure
subplot (2 ,1 ,1)
plot (t , bit_mod );
xlabel ('Time')
ylabel ('Bit Value')
title ('Original Bit Stream')
xlim ([0 2])
ylim ([ -0.1 1.1])

subplot (2 ,1 ,2)
stem ( t (1: Tb / dt : length ( t ) -1) , bpsk_sym ) ;
xlabel ('Time')
ylabel ('Amplitude')
title ('BPSK Symbol Train')
xlim ([0 2])
ylim ([ -1.1 1.1])

% plotting the BPSK symbol train
figure
stem ( bpsk_sym ) ;
xlabel ('n')
ylabel ('Bit Value')
title ('BPSK Symbol Train')
xlim ([0 50])
ylim ([ -1.5 1.5])

%% Convolving with the sinc /raised - cosine function




alpha_1 = 0.5;
alpha_2 = 1;
ts = -L*Tb :dt: L*Tb;

Rx_sinc = sinc(fs*ts);
Rx_rcos1 = rcos(alpha_1,fs,ts);
Rx_rcos2 = rcos(alpha_2,fs,ts);

tx_stream_Sinc = conv(BPSK_time ,Rx_sinc ,'same') ;
tx_stream_Rcos1 = conv(BPSK_time ,Rx_rcos1 ,'same') ;
tx_stream_Rcos2 = conv(BPSK_time ,Rx_rcos2 ,'same') ;

figure
plot (t , tx_stream_Sinc ) ;
hold on
stem ( t (1: Tb / dt : length ( t ) -1) , bpsk_sym ,'k','LineStyle','--') ;
xlabel ('Time')
ylabel ('Amplitude')
title ('Received Signal - Sinc ( Without Noise )')
xlim ([0 1])
grid on

figure
plot (t , tx_stream_Rcos1 ) ;
hold on
stem ( t (1: Tb / dt : length ( t ) -1) , bpsk_sym ,'k','LineStyle','--') ;
xlabel ('Time')
ylabel ('Amplitude')
title ('Received Signal - Raised Cosine with \ alpha = 0.5 ( Without Noise )')
xlim ([0 1])
grid on

figure
plot (t , tx_stream_Rcos2 ) ;
hold on
stem ( t (1: Tb / dt : length ( t ) -1) , bpsk_sym ,'k','LineStyle','--') ;
xlabel ('Time')
ylabel ('Amplitude')
title ('Received Signal - Raised Cosine with \ alpha = 1 ( Without Noise )')
xlim ([0 1])
grid on

%% Eye Diagrams

eyediagram ( tx_stream_Sinc ,2*(1/ fs ) / dt ,2/ fs ) ;
eyediagram ( tx_stream_Rcos1 ,2*(1/ fs ) / dt ,2/ fs ) ;
eyediagram ( tx_stream_Rcos2 ,2*(1/ fs ) / dt ,2/ fs ) ;


%% Adding Noise

Eb = 1; % bit energy
gamma = 10; % (Eb/N0) in dB
N0 = Eb *10^( -0.1* gamma ) ; % noise power spectral density

awgn = sqrt ( N0 /2) * randn (1 , length ( tx_stream_Sinc ) ) ;

rx_noise_Sinc = tx_stream_Sinc + awgn ;
rx_noise_Rcos1 = tx_stream_Rcos1 + awgn ;
rx_noise_Rcos2 = tx_stream_Rcos2 + awgn ;

% plotting the received signals with noise
figure
plot (t , rx_noise_Sinc ) ;
hold on
stem ( t (1: Tb / dt : length ( t ) -1) , bpsk_sym ,'k','LineStyle','--') ;
xlabel ('Time (s)')
ylabel ('Amplitude ')
title ( strcat ('Received Signal (E_b/N_0 = ',num2str ( gamma ) ,'B) - Sinc ') )
xlim ([0 1])
grid on

figure
plot (t , rx_noise_Rcos1 ) ;
hold on
stem ( t (1: Tb / dt : length ( t ) -1) , bpsk_sym ,'k','LineStyle','--') ;
xlabel ('Time (s)')
ylabel ('Amplitude ')
title ( strcat ('Received Signal (E_b/N_0 = ',num2str ( gamma ) ) )
xlim ([0 1])
grid on

figure
plot (t , rx_noise_Rcos2 ) ;
hold on
stem ( t (1: Tb / dt : length ( t ) -1) , bpsk_sym ,'k','LineStyle','--') ;
xlabel ('Time (s)')
ylabel ('Amplitude ')
title ( strcat ('Received Signal (E_b/N_0 = ',num2str ( gamma )  ))
xlim ([0 1])
grid on

%% Eye Diagrams with Noise

eyediagram ( rx_noise_Sinc ,2*(1/ fs ) / dt ,2/ fs ) ;
eyediagram ( rx_noise_Rcos1 ,2*(1/ fs ) / dt ,2/ fs ) ;
eyediagram ( rx_noise_Rcos2 ,2*(1/ fs ) / dt ,2/ fs ) ;



clear ; close all;
L  = 1000000;
fs = 10;            % Sampling Frequency
Tb = 1/fs;
t =0:Tb:(L-1)*Tb;

% CREATING RANDOM BINARY SEQUENCE
data = (rand(1,L)>0.5);

figure
stem(data,'LineStyle',"--");
grid("on")
title("Random Binary Sequence")
ylabel("bit-value")
ylim([-0.5 1.5]);
xlim([0 20]);

% Map the Sequence Using BPSK
x = 2*data-1;

figure;
stem(t,x,'LineStyle',"--");
grid("on")
title("Transmitted Pulse")
ylabel("Amplitude")
xlabel("Time(s)")
ylim([-1.2 1.2]);
xlim([0 2]);

% Channel-Delay Characteristics
h = [0.3 0.9 0.4];
l_ = floor(length(h)/2);                % Number of samples to disgard
y = conv(x,h,"same");


figure
hold on;
stem(t,x,'LineStyle',"none");
stem(t,y,'filled',"LineStyle","--");
grid("on")
title("Received Signal (noise-free)")
ylabel("Amplitude")
xlabel("Time(s)")
ylim([-3 3]);
xlim([0 2]);
legend(["Transmitted","Received"])

% Channel Noise
noisy_y = awgn(y,3.14);

figure
hold on;
stem(t,x,'LineStyle',"none");
stem(t,noisy_y,'filled',"LineStyle","-.");
grid("on")
title("Received Signal")
ylabel("Amplitude")
xlabel("Time(s)")
ylim([-3 3]);
xlim([0 2]);
legend(["Transmitted","Received"])


% M-Tap Equalization Filter
figure
stem(t,x,'LineStyle',"none");
hold on;

for M = [3 9]
    he = ZF_coff(max(M,length(h)),h,L);
    recovered = conv(y,he([(end-M)/2+1:(end+M)/2]),"same");
    
    stem(t,recovered,'filled',"LineStyle","--");
end

grid("on")
title("Recovered Signal (Noise-Free)")
ylabel("Amplitude")
xlabel("Time(s)")
ylim([-2 2]);
xlim([0 2]);
legend(["Transmitted","Recovered M = 3","Recovered M=9"])

% Equalization on Noisy Signal
recovered = conv(noisy_y,he([(end-M)/2+1:(end+M)/2]),"same");


figure
stem(t,x,'LineStyle',"none");
hold on;
stem(t,recovered,'filled',"LineStyle","-.");

grid("on")
title("Recovered Signal")
ylabel("Amplitude")
xlabel("Time(s)")
ylim([-3 3]);
xlim([0 2]);
legend(["Transmitted","Recovered"])


% De-Modulation
d_mod = recovered > 0;

figure
hold on;
stem(t,data,'LineStyle',"none");
stem(t,d_mod,'filled',"LineStyle","-.");
grid("on")
title("Demodulated Data")
ylabel("Value")
xlabel("Time(s)")
ylim([-0.5 1.2]);
xlim([0 2]);
legend(["Original","Demodulated"],'Location','southeast')

% BER calculation
n_errors = biterr(data,d_mod);
p_error = n_errors/L

figure;
hold on;
for snr = [0 1 2 3 4 5 6 7 8 9 10]+3.14
    noisy_y = awgn(y,snr);
    p_error = [];
    for M = 1:2:9
        he = ZF_coff(max(M,length(h)),h,L);
        recovered = conv(noisy_y,he([(end-M)/2+1:(end+M)/2]),"same");
        d_mod = recovered > 0;
        n_errors = biterr(data,d_mod);
        p_error = [p_error n_errors/L];
    end
    plot(1:2:9,p_error);
end
grid on;
title("Error Performance of M-Tapped Filter");
ylabel("BER")
xlabel("M")
legend(["E_b/N_0 = 0dB","E_b/N_0 = 1dB","E_b/N_0 = 2dB","E_b/N_0 = 3dB","E_b/N_0 = 4dB","E_b/N_0 = 5dB","E_b/N_0 = 6dB","E_b/N_0 = 7dB","E_b/N_0 = 8dB","E_b/N_0 = 9dB","E_b/N_0 = 10dB"],"location",'northwest')

Eb_N0 = 9:.5:15;
figure;
for M = 3:2:9
    p_error = [];
    for snr = Eb_N0 + 3.14
        noisy_y = awgn(y,snr); 
        he = ZF_coff(max(M,length(h)),h,L);
        recovered = conv(noisy_y,he([(end-M)/2+1:(end+M)/2]),"same");
        d_mod = recovered > 0;
        n_errors = biterr(data,d_mod);
        p_error = [p_error n_errors/L];
    end
    semilogy(Eb_N0,p_error)
    hold on;    
end
grid on;
axis([9 15 0 0.1])
yticks(0:0.04:0.1);
title("Error Performance - SNR");
xlabel("{E_b}/{N_0}(dB Scale)");
ylabel("BER(dB Scale)");
legend(["M=3","M=5","M=7","M=9"]);

function [output] = rcos(alpha,Rb,t)
    output = cos(pi*alpha*Rb*t).*sinc(Rb*t)./(1-(2*alpha*Rb*t).^2);
end

function [he] = ZF_coff(M,hc,L) 
    hc = [zeros(1,(M-length(hc))/2) hc zeros(1,(M-length(hc))/2)];
    delta = zeros(1,M);
    delta(round(end/2)) = 1;
    p = conv(delta,hc);
    r = p([ceil(end/2):-1:1]);
    c = p([ ceil(end/2):end ]);
    T = toeplitz(c,r);
    coff = inv(T)*reshape(delta,[M 1]);
    he = reshape(coff,[1,M]);
end