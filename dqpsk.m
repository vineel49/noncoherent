% Noncoherent detection of DQPSK in AWGN

clear all
close all
clc
SNR_dB = 12; % SNR per bit in dB
num_bit = 10^5;
tic()
% source
a = randi([0 1],1,num_bit);

% Differential encoding
dqpsk_seq = zeros(1,num_bit/2);
ref_sym = 1+1i; % reference symbol
dqpsk_rules = [exp(1i*0) exp(1i*pi/2) exp(1i*3*pi/2) exp(1i*pi)]; % dqpsk encoding rules
dqpsk_seq(1) = ref_sym * dqpsk_rules(2*a(1)+a(2)+1);
for i = 2:(num_bit/2) 
   dqpsk_seq(i) = dqpsk_seq(i-1)* dqpsk_rules( 2*a(2*i-1)+a(2*i)+1);
end

% AWGN
SNR = 10^(0.1*SNR_dB);
noise_var_1D = 0.5*2/(2*SNR);
noise =  normrnd(0,sqrt(noise_var_1D),1,num_bit/2)+1i*normrnd(0,sqrt(noise_var_1D),1,num_bit/2);

% Channel output
Chan_Op = dqpsk_seq + noise;

% Noncoherent Receiver
comp_mat = zeros(4,num_bit/2-1);
comp_mat(1,1:2:end) = real( conj(Chan_Op(1:2:end-1)).* Chan_Op(2:2:end)*exp(-1i*0));
comp_mat(2,1:2:end) = real( conj(Chan_Op(1:2:end-1)).* Chan_Op(2:2:end)*exp(-1i*pi/2));
comp_mat(3,1:2:end) = real( conj(Chan_Op(1:2:end-1)).* Chan_Op(2:2:end)*exp(-1i*3*pi/2));
comp_mat(4,1:2:end) = real( conj(Chan_Op(1:2:end-1)).* Chan_Op(2:2:end)*exp(-1i*pi));   
    
comp_mat(1,2:2:end) = real( conj(Chan_Op(2:2:end-1)).* Chan_Op(3:2:end)*exp(-1i*0));
comp_mat(2,2:2:end) = real( conj(Chan_Op(2:2:end-1)).* Chan_Op(3:2:end)*exp(-1i*pi/2));
comp_mat(3,2:2:end) = real( conj(Chan_Op(2:2:end-1)).* Chan_Op(3:2:end)*exp(-1i*3*pi/2));
comp_mat(4,2:2:end) = real( conj(Chan_Op(2:2:end-1)).* Chan_Op(3:2:end)*exp(-1i*pi));

% comparing column wise
[~,indices] = max(comp_mat,[],1);
% MAPPING INDICES TO QPSK SYMBOLS
QPSK_SYM = [1+1i 1-1i -1+1i -1-1i];
temp = QPSK_SYM(indices);
% DEMAPPING QPSK SYMBOLS TO BITS
dec_a = zeros(1,num_bit-2);
dec_a(1:2:end) = real(temp)<0;
dec_a(2:2:end) = imag(temp)<0;

% Bit error rate
BER = nnz(a(3:end)- dec_a)/(num_bit-2)
toc()