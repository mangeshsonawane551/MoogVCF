classdef moogVCF_TrapezoidalSounds_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        rResonantFrequencyKnobLabel  matlab.ui.control.Label
        rResonantFrequencyKnob       matlab.ui.control.Knob
        f0FundamentalFrequencyKnobLabel  matlab.ui.control.Label
        f0FundamentalFrequencyKnob   matlab.ui.control.Knob
        stopButton                   matlab.ui.control.Button
        PlaySoundButton              matlab.ui.control.Button
        BasicMOOGVCFLabel            matlab.ui.control.Label
    end

    methods (Access = private)

        % Value changing function: rResonantFrequencyKnob
        function rResonantFrequencyKnobValueChanging(app, event)
            changingValue = event.Value;
            app.rResonantFrequencyKnob.Value=changingValue;
        end

        % Value changing function: f0FundamentalFrequencyKnob
        function f0FundamentalFrequencyKnobValueChanging(app, event)
      
            changingValue = event.Value;
            app.f0FundamentalFrequencyKnob.Value=changingValue;
            
  

        end

        % Button pushed function: stopButton
        function stopButtonPushed(app, event)
          
clear all;
        end

        % Button pushed function: PlaySoundButton
        function PlaySoundButtonPushed(app, event)
            
                        Fs = 44100; 

% Time step
k = 1/Fs; 

% time step
f0 = app.f0FundamentalFrequencyKnob.Value;

% Resonant frequency in Hz
r = app.rResonantFrequencyKnob.Value; 
% Error checking for r
if r > 1 || r < 0
   error('r must be between 0 and 1 only');
end

% Angular Frequency 
w0 = 2*pi*f0;

% Simulation duration in sec
Tf = 2;

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

%tt= 0:Nf-1;
% response vector
t = 0:k:Tf-k;
un = sawtooth(2*pi*f0*t);

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

soundsc(yTrapOut,Fs);
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
            app.rResonantFrequencyKnobLabel.FontWeight = 'bold';
            app.rResonantFrequencyKnobLabel.Position = [101 169 131 30];
            app.rResonantFrequencyKnobLabel.Text = {'r'; '(Resonant Frequency)'};

            % Create rResonantFrequencyKnob
            app.rResonantFrequencyKnob = uiknob(app.UIFigure, 'continuous');
            app.rResonantFrequencyKnob.Limits = [0.1 0.9];
            app.rResonantFrequencyKnob.ValueChangingFcn = createCallbackFcn(app, @rResonantFrequencyKnobValueChanging, true);
            app.rResonantFrequencyKnob.Position = [99 233 136 136];
            app.rResonantFrequencyKnob.Value = 0.3;

            % Create f0FundamentalFrequencyKnobLabel
            app.f0FundamentalFrequencyKnobLabel = uilabel(app.UIFigure);
            app.f0FundamentalFrequencyKnobLabel.HorizontalAlignment = 'center';
            app.f0FundamentalFrequencyKnobLabel.Position = [395.5 177 142 28];
            app.f0FundamentalFrequencyKnobLabel.Text = {'f0'; '(Fundamental Frequency)'};

            % Create f0FundamentalFrequencyKnob
            app.f0FundamentalFrequencyKnob = uiknob(app.UIFigure, 'continuous');
            app.f0FundamentalFrequencyKnob.Limits = [50 20000];
            app.f0FundamentalFrequencyKnob.ValueChangingFcn = createCallbackFcn(app, @f0FundamentalFrequencyKnobValueChanging, true);
            app.f0FundamentalFrequencyKnob.Position = [403 239 127 127];
            app.f0FundamentalFrequencyKnob.Value = 1000;

            % Create stopButton
            app.stopButton = uibutton(app.UIFigure, 'push');
            app.stopButton.ButtonPushedFcn = createCallbackFcn(app, @stopButtonPushed, true);
            app.stopButton.IconAlignment = 'center';
            app.stopButton.BackgroundColor = [1 1 1];
            app.stopButton.FontWeight = 'bold';
            app.stopButton.Position = [417 93 100 23];
            app.stopButton.Text = 'stop';

            % Create PlaySoundButton
            app.PlaySoundButton = uibutton(app.UIFigure, 'push');
            app.PlaySoundButton.ButtonPushedFcn = createCallbackFcn(app, @PlaySoundButtonPushed, true);
            app.PlaySoundButton.BackgroundColor = [1 1 1];
            app.PlaySoundButton.FontWeight = 'bold';
            app.PlaySoundButton.Position = [117 93 100 23];
            app.PlaySoundButton.Text = ' Play  Sound';

            % Create BasicMOOGVCFLabel
            app.BasicMOOGVCFLabel = uilabel(app.UIFigure);
            app.BasicMOOGVCFLabel.HorizontalAlignment = 'center';
            app.BasicMOOGVCFLabel.FontName = 'Chalkboard SE';
            app.BasicMOOGVCFLabel.FontSize = 33;
            app.BasicMOOGVCFLabel.FontWeight = 'bold';
            app.BasicMOOGVCFLabel.Position = [177.5 403 274 51];
            app.BasicMOOGVCFLabel.Text = 'Basic MOOG VCF';
        end
    end

    methods (Access = public)

        % Construct app
        function app = S1889125_moogVCF_TrapezoidalSounds_GUI

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
