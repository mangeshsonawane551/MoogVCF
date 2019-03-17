%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Author: MANGESH SONAWANE
% Program details:It is a program that simulates the impulse response of a 
% Moog VCF ladder filter based on forward and backward integrators of the 
% linear state space system and plots the transfer function of the system,
% forward and backward output signal transfer function
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clear all;
clc; 
tic
% sample rate
Fs = 44100;
if isnumeric(Fs)==0
    error('Sample rate should be numeric');
end

% Time step
k = 1/Fs; 

% Resonant frequency in Hz
f0 = 1000; 
if f0<20 || f0>20000
    error('Please enter fundamental within 20-20Khz range');
end

% Angular Frequency 
w0 = 2*pi*f0; 

% tuning parameter 
r = 0.3;
% Error checking for r
if r > 1 || r < 0
   error('r must be between 0 and 1 only');
end

% Simulation duration
Tf = 3;

% Total number of frames/samples
Nf = floor(Tf*Fs);

% System matrix 'A'
A = w0*[-1 0 0 -4*r; 1 -1 0 0; 0 1 -1 0; 0 0 1 -1]; 
if real(eig(A))>0
    error('Eigen values are not in the range')
end

% Forcing vector 'b'
b = w0*[1; 0; 0; 0];

% State mixture vector 'c'
C = transpose([0; 0; 0; 1]); 

% Unit impulse response vector
un = zeros(Nf,1); 
un(1) = 1;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%            Implementation of forward integration method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initialize state vectors

% Current time step
xf2 = zeros(4,1); 
% Current-1 time step
xf1 = zeros(4,1); 

% Initialize output vector
yForward = zeros(Nf,1);

% Set coefficient
I = eye(4);
coefficient = I + k*A;

% Calculate stability condition
 if k > -2/(w0*(-1 + sqrt(2)*r^(1/4)*exp(1i*(pi/4 + pi/2))))
           error('|1 + k*eig (A) | ? 1 ' );
 end

%Forward equation implementation
for i = 2:Nf
   
   xf2 = coefficient*xf1 + k*b*un(i-1);
   yForward(i) = C*xf2;
   
   % update value
   xf1 = xf2;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Implementation of  backward integration method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initializing state vectors

% Current time step
xb2 = zeros(4,1); 

% Current-1 time step
xb1 = zeros(4,1); 

% Initialize output vector
yBackward = zeros(Nf,1);

% Compute coefficient Q
Q = eye(4) - k*A;

% Backward equation implementation
for i = 1:Nf
  
   d = xb1 + k*b*un(i);
   xb2 = Q\d;
   yBackward(i) = C*xb2;
   
   % Update values
   xb1 = xb2;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Exact transfer function of system
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Maxmimun frequency
fNyquist = Fs/2;

% Number of bin with 1000 spacing
fbin = 0:fNyquist/1000:fNyquist;

%Length of bins
flen = length(fbin);

%Initialise Hexact vector
Hexact = zeros(flen,1);

% Calculating Hexact of system
for i = 1:flen
    x2exact = (2*pi*1i*fbin(i)*I-A)\b;
    Hexact(i) = C*x2exact;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Transfer functions of output
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Calculating Fourier transforms

%DFT of input signal
UNFFT = fft(un);

% DFT of output signal by forward method 
YForwardFFT = fft(yForward);

% DFT of output signal by backward method
YBackwardFFT = fft(yBackward);


% Transfer function for forward method
HForwardT = UNFFT./YForwardFFT;

% Transfer function for backward method
HBackwardT = UNFFT./YBackwardFFT;

% log conversion 
HForwardLog = log(HForwardT) - max(log(HForwardT));
HBackwardLog = log(HBackwardT) - max(log(HBackwardT));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                                   %PLOT
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Sample length
n = 0:Nf-1; 

% frequency bins [Hz]
fk = n'*Fs/Nf; 
figure(1);
plot(fk,abs(HBackwardLog),'r','Linewidth', 1);
hold on
plot(fk,abs(HForwardLog), 'm','Linewidth', 1);
hold on
plot(fbin,real(log(Hexact)), 'b','Linewidth', 1);
xlim([0 fNyquist]);
title('Transfer Functions');
xlabel('Frequency (Hz)'); 
ylabel('Magnitude');
legend('Forward', 'Backward', 'Exact')
grid on;

toc
