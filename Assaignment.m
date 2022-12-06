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

