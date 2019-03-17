classdef moogVCF_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        rResonantFrequencyKnobLabel    matlab.ui.control.Label
        rResonantFrequencyKnob         matlab.ui.control.Knob
        f0102FundamentalfrequencyKnobLabel  matlab.ui.control.Label
        f0102FundamentalfrequencyKnob  matlab.ui.control.Knob
        UIAxes                         matlab.ui.control.UIAxes
        HexactButton                   matlab.ui.control.Button
        HForwardButton                 matlab.ui.control.Button
        HbackwardButton                matlab.ui.control.Button
        HTrapezoidButton               matlab.ui.control.Button
        PlotLabel                      matlab.ui.control.Label
        TransferFunctionforMoogVCFLabel  matlab.ui.control.Label
    end

    methods (Access = private)

        % Value changed function: rResonantFrequencyKnob
        function rResonantFrequencyKnobValueChanged(app, event)
            value = app.rResonantFrequencyKnob.Value;
            
        end

        % Value changed function: f0102FundamentalfrequencyKnob
        function f0102FundamentalfrequencyKnobValueChanged(app, event)
            value = app.f0102FundamentalfrequencyKnob.Value;
            
        end

        % Button pushed function: HexactButton
        function HexactButtonPushed(app, event)
    Fs = 44100;

% Time step
k = 1/Fs; 

% Resonant frequency in Hz
f0 = (app.f0102FundamentalfrequencyKnob.Value)*10^2; 

% Angular Frequency 
w0 = 2*pi*f0; 

% tuning parameter 
r = app.rResonantFrequencyKnob.Value;
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
% Maxmimun frequency
fNyquist = Fs/2;

I = eye(4);
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

plot(app.UIAxes,fbin,real(log(Hexact)), 'b','LineWidth',1);
xlim(app.UIAxes,[0 fNyquist]);
legend(app.UIAxes,'Hexact');

        end

        % Button pushed function: HbackwardButton
        function HbackwardButtonPushed(app, event)
            Fs = 44100;

% Time step
k = 1/Fs; 

% Resonant frequency in Hz
f0 = (app.f0102FundamentalfrequencyKnob.Value)*10^2;

% Angular Frequency 
w0 = 2*pi*f0; 

% tuning parameter 
r = app.rResonantFrequencyKnob.Value;
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

%DFT of input signal
UNFFT = fft(un);


% DFT of output signal by backward method
YBackwardFFT = fft(yBackward);



% Transfer function for backward method
HBackwardT = UNFFT./YBackwardFFT;


HBackwardLog = log(HBackwardT) - max(log(HBackwardT));

% Sample length
n = 0:Nf-1; 

% frequency bins [Hz]
fk = n'*Fs/Nf; 
fNyquist = Fs/2;


plot(app.UIAxes, fk,abs(HBackwardLog),'r','Linewidth', 1);
xlim(app.UIAxes,[0 fNyquist]);
legend(app.UIAxes,'HBackward');
        end

        % Button pushed function: HTrapezoidButton
        function HTrapezoidButtonPushed(app, event)
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% sample rate
Fs = 44100; 

% Time step
k = 1/Fs; 

% time step
f0 = (app.f0102FundamentalfrequencyKnob.Value)*10^2;

% Resonant frequency in Hz
r = app.rResonantFrequencyKnob.Value; 
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


% Maxmimun frequency
fnyquist = Fs/2;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                                   Plot 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Sample length
n = 0:Nf-1; 



% Frequency bins [Hz]
fk = n'*Fs/Nf; 


plot(app.UIAxes, fk, abs(HTrapLog),'m','LineWidth',1);

xlim(app.UIAxes, [0 fnyquist]);

legend(app.UIAxes, 'Trapezoidal')



toc 


        end

        % Button pushed function: HForwardButton
        function HForwardButtonPushed(app, event)
            % sample rate
Fs = 44100;

% Time step
k = 1/Fs; 

% Resonant frequency in Hz
f0 = 1000; 

% Angular Frequency 
w0 = 2*pi*f0; 

% tuning parameter 
r = app.rResonantFrequencyKnob.Value;
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

%DFT of input signal
UNFFT = fft(un);

% DFT of output signal by forward method 
YForwardFFT = fft(yForward);




% Transfer function for forward method
HForwardT = UNFFT./YForwardFFT;


fNyquist = Fs/2;
% log conversion 
HForwardLog = log(HForwardT) - max(log(HForwardT));
n = 0:Nf-1; 
fk = n'*Fs/Nf; 
plot(app.UIAxes,fk,abs(HForwardLog), 'm','LineWidth',1);
xlim(app.UIAxes,[0 fNyquist]);
legend(app.UIAxes,'HForward');
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'UI Figure';

            % Create rResonantFrequencyKnobLabel
            app.rResonantFrequencyKnobLabel = uilabel(app.UIFigure);
            app.rResonantFrequencyKnobLabel.HorizontalAlignment = 'center';
            app.rResonantFrequencyKnobLabel.FontSize = 11;
            app.rResonantFrequencyKnobLabel.FontWeight = 'bold';
            app.rResonantFrequencyKnobLabel.Position = [25 264 120 26];
            app.rResonantFrequencyKnobLabel.Text = {'r '; '(Resonant Frequency)'};

            % Create rResonantFrequencyKnob
            app.rResonantFrequencyKnob = uiknob(app.UIFigure, 'continuous');
            app.rResonantFrequencyKnob.Limits = [0.1 0.9];
            app.rResonantFrequencyKnob.ValueChangedFcn = createCallbackFcn(app, @rResonantFrequencyKnobValueChanged, true);
            app.rResonantFrequencyKnob.Position = [55 313 60 60];
            app.rResonantFrequencyKnob.Value = 0.5;

            % Create f0102FundamentalfrequencyKnobLabel
            app.f0102FundamentalfrequencyKnobLabel = uilabel(app.UIFigure);
            app.f0102FundamentalfrequencyKnobLabel.HorizontalAlignment = 'center';
            app.f0102FundamentalfrequencyKnobLabel.FontSize = 11;
            app.f0102FundamentalfrequencyKnobLabel.FontWeight = 'bold';
            app.f0102FundamentalfrequencyKnobLabel.Position = [30 88 136 26];
            app.f0102FundamentalfrequencyKnobLabel.Text = {'f0 (10^2)'; '(Fundamental frequency)'};

            % Create f0102FundamentalfrequencyKnob
            app.f0102FundamentalfrequencyKnob = uiknob(app.UIFigure, 'continuous');
            app.f0102FundamentalfrequencyKnob.Limits = [10 200];
            app.f0102FundamentalfrequencyKnob.ValueChangedFcn = createCallbackFcn(app, @f0102FundamentalfrequencyKnobValueChanged, true);
            app.f0102FundamentalfrequencyKnob.Position = [66 148 60 60];
            app.f0102FundamentalfrequencyKnob.Value = 50;

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Transfer function Plot')
            xlabel(app.UIAxes, 'Frequency')
            ylabel(app.UIAxes, 'Magnitude')
            app.UIAxes.FontWeight = 'bold';
            app.UIAxes.Position = [169 56 472 363];

            % Create HexactButton
            app.HexactButton = uibutton(app.UIFigure, 'push');
            app.HexactButton.ButtonPushedFcn = createCallbackFcn(app, @HexactButtonPushed, true);
            app.HexactButton.Position = [82 32 100 22];
            app.HexactButton.Text = 'H exact';

            % Create HForwardButton
            app.HForwardButton = uibutton(app.UIFigure, 'push');
            app.HForwardButton.ButtonPushedFcn = createCallbackFcn(app, @HForwardButtonPushed, true);
            app.HForwardButton.Position = [191 32 100 22];
            app.HForwardButton.Text = 'H Forward';

            % Create HbackwardButton
            app.HbackwardButton = uibutton(app.UIFigure, 'push');
            app.HbackwardButton.ButtonPushedFcn = createCallbackFcn(app, @HbackwardButtonPushed, true);
            app.HbackwardButton.Position = [300 32 100 22];
            app.HbackwardButton.Text = 'H backward';

            % Create HTrapezoidButton
            app.HTrapezoidButton = uibutton(app.UIFigure, 'push');
            app.HTrapezoidButton.ButtonPushedFcn = createCallbackFcn(app, @HTrapezoidButtonPushed, true);
            app.HTrapezoidButton.Position = [410 32 100 22];
            app.HTrapezoidButton.Text = 'H Trapezoid';

            % Create PlotLabel
            app.PlotLabel = uilabel(app.UIFigure);
            app.PlotLabel.FontSize = 18;
            app.PlotLabel.FontWeight = 'bold';
            app.PlotLabel.Position = [41 29 42 25];
            app.PlotLabel.Text = 'Plot';

            % Create TransferFunctionforMoogVCFLabel
            app.TransferFunctionforMoogVCFLabel = uilabel(app.UIFigure);
            app.TransferFunctionforMoogVCFLabel.HorizontalAlignment = 'center';
            app.TransferFunctionforMoogVCFLabel.FontName = 'Athelas';
            app.TransferFunctionforMoogVCFLabel.FontSize = 15;
            app.TransferFunctionforMoogVCFLabel.FontWeight = 'bold';
            app.TransferFunctionforMoogVCFLabel.Position = [-98 418 837 45];
            app.TransferFunctionforMoogVCFLabel.Text = 'Transfer Function for Moog VCF using forward backward and trapezoidal integration methods';
        end
    end

    methods (Access = public)

        % Construct app
        function app = S1889125_moogVCF_GUI

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
