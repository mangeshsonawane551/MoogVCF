%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Author: MANGESH SONAWANE
% Program details: It is a program that simulates the impulse response of a 
% Moog VCF ladder filter based on trapezoidal integration of the 
% linear state space system and plots the transfer function of the system,
% and trapezoidal output signal transfer function
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
tic;
clear all;
clc;

% sample rate
Fs = 44100; 
if isnumeric(Fs)==0
    error('Sample rate should be numeric');
end
% Time step
k = 1/Fs; 

% Fundamental frequency
f0 = 1000; 
if f0<20 || f0>20000
    error('Please enter fundamental within 20-20Khz range');
end

% Resonant frequency in Hz
r = 0.5; 
% Error checking for r
if r > 1 || r < 0
   error('r must be between 0 and 1 only');
end

% Angular Frequency 
w0 = 2*pi*f0;

% Simulation duration in sec
Tf = 3;

% Total number of frames/ samples
Nf = floor(Tf*Fs); 

% System matrix 'A'
A = w0*[-1 0 0 -4*r; 1 -1 0 0; 0 1 -1 0; 0 0 1 -1]; 
if real(eig(A))>0
    error('Eigen values are not in the range')
end

% Forcing vector 'b'
b = w0*[1; 0; 0; 0]; 

% State mixture vector 'c'
c = transpose([0; 0; 0; 1]);

% Unit impulse response vector
un = zeros(Nf,1); 
un(1) = 1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%            Implementation of trapezoidal integration method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initialize state vectors

% Current time step
xT2 = zeros(4,1); 
% Current-1 time step
xT1 = zeros(4,1); 

% Initialize output vector
yTrapOut = zeros(Nf,1);

% Set coefficient
coefficient1 = eye(4) + k*A/2;
coefficient2 = eye(4) - k*A/2;

% Trapezoidal integration equation implementation
for n = 1:Nf
   
   
   d = coefficient1*xT1 +k*coefficient1*b*un(n);
   xT2 = coefficient2\d;
   yTrapOut(n) = c*xT2;
   
   % Update value
   xT1 = xT2;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Transfer functions of output
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Calculating Fourier transforms

%DFT of input signal
UNFFT = fft(un);

% DFT of output signal 
YTrapFFT = fft(yTrapOut);

% Transfer functions
HTrap = UNFFT./YTrapFFT;

%Log Conversion
HTrapLog = log(HTrap) - max(log(HTrap));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                   Exact transfer function of system
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Maxmimun frequency
fnyquist = Fs/2;

% Number of bin with 1000 spacing
fbin = 0:fnyquist/1000:fnyquist;

%Length of bins
flen = length(fbin);

%Initialise Hexact vector
Hexact = zeros(flen,1);


% Calculating Hexact of system
for n = 1:flen
    xexact = (2*pi*1i*fbin(n)*eye(4)-A)\b;
    Hexact(n) = c*xexact;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                                   Plot 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Sample length
n = 0:Nf-1; 

% Frequency bins [Hz]
fk = n'*Fs/Nf; 

figure(1);
plot(fk, abs(HTrapLog),'r','Linewidth', 1);
hold on
plot(fbin, real(log(Hexact)), 'm','Linewidth', 1);
xlim([0 fnyquist]);
title('Transfer Functions Moog vcf Trapezoidal integration');
xlabel('Frequency (Hz)'); 
ylabel('Magnitude');
legend('Trapezoidal', 'Exact')
grid on;


toc 

