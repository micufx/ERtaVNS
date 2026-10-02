%% runDAQMethodsOfLimits
  

% time is set to 1 minute 
% zero_signal removed so we dont have to wait 13s for the stimulus 
% base script for NiDAQ-Analog Data Generation
%-Heiko Stecher 2018

%% tidy up
clear all 
close all
clc

%% NiDAQ parametres
deviceID          = 'Dev2';          % device ID of the nidaq, check in NiMax, different on every computer, NiMax opens automatically if device is plugged in
Fs                = 200000;             % sampling rate, with 500k a precise pulsewidth is created

%% NOTE TO FS 

%DO NOT USE THE FULL 500kHz. The signal is unnecessary precise and too
%large. Likely your PC will not have enough memory to generate a full stim
%sequence. The method below I used is very crude and has 20micro seconds
%variance. This can easily be improved without increasing the sampling rate
%of the carrier frequency. The problem is likely undersampling at 10kHz.
%With a bit of math you can figure out which frequency works best for your
%given duty cycle. 

%% Signal parameters
Stim_Freq = 25; % Stimulation with 25Hz Pulses 
%Stim_window_len = 3.4; % we chose a stimulation window of 3.4 seconds
Stim_window_len = 5; % we chose a stimulation window of 3.4 seconds
t_stim = 0:1/Fs:1*Stim_window_len - (1/Fs); % calculate the time vector for the stimulation sequence
%t_stim_dur = 1; % total stim´duration in minutes we try x minutes
%t_stim_dur_s = t_stim_dur*60; % recalculate to seconds
t_stim_dur_s = 3; % recalculate to seconds
p_width = 250e-6; %250 microseconds pulsewidth 
%% Amplitude input
amp = input('Type in amplitude: ', 's');  % input as string

if contains(amp, ',')
    error('Use a dot (.) as the decimal separator, not a comma (,).');
end
if isnan(str2double(amp))
    error('You did not specify an amplitude!');
else
    amp = str2double(amp);
end

if amp > 5
    error('The amplitude is too big! The maximum amplitude is 5 mA');
end

%% find out how often we have to repeat the signal to fit it to the total stimulation duration
sequence_time = Stim_window_len; %length of the window to stimulate % duration in seconds of one signal cycle (used for ramp)
stim_reps = round(t_stim_dur_s/sequence_time); %number f times we repeat the signal

%% PRÄZISERE ART DER STIMULUS GENERIERUNG
cycle_time = (1/Stim_Freq);
t = -cycle_time/2:1/Fs:cycle_time/2-(1/Fs); %von minus eine halbe Periode bis plus eine halbe Periode
w = p_width; %Pulsbreite
x = amp*rectpuls(t,w); % Generiere einen Pulse
stim_unipol2 = repmat(x,1,(round(Stim_window_len/(1/Stim_Freq)))); %repmat pulse % generates array with x as value of the puls, 1 row number, rounds numbers of repeats that fit in 3.4s with 25 Hz
% plot(x)


duty_mono = round(dutycycle(stim_unipol2,t_stim)*100,2);
fprintf('Minimum duty cycle was %f %%, mean duty cycle was %f %%, and max duty cycle was %f %%. \n',min(duty_mono), mean(duty_mono), max(duty_mono))
% calculate the pulse width for all signals
pulse_width_mono = pulsewidth(stim_unipol2,t_stim)*1000000;
fprintf('Minimum pulse width was %f %%, mean pulse width was %f %%, and max pulse width was %.f %%. \n',min(pulse_width_mono), mean(pulse_width_mono), max(pulse_width_mono))

signal = stim_unipol2;

% create a signal of a given time (min) with the fitting number of
% stim_reps

signal_min = repmat(signal,1,stim_reps )'; % generates a signal of time t_stim_dur(60minutes as an example) 

% plot to see the result 
% t_plot = 0:1/Fs:(stim_reps*sequence_time-1/Fs);
% plot(t_plot, signal_min (1:length(t_plot)));
% axis([0 (stim_reps*sequence_time-1/Fs) -10 10]);
% pause(1); % show figure for 5 seconds
% close; % close figure

%% initializing nidaq
daqreset;                                 % In case there was still an open session, reset the NiDAQ
s = daq.createSession('ni');              % open session % https://de.mathworks.com/help/daq/ref/daq.createsession.html
s.addAnalogOutputChannel(deviceID,[0], 'Voltage'); % Opens one output channel. Note: first channel is designated "0" %https://de.mathworks.com/help/daq/ref/addanalogoutputchannel.html
set(s, 'Rate', Fs);                       % Set sampling rate
s.IsContinuous = true;                    % Sets NiDAQ to run until stopped % https://de.mathworks.com/help/daq/ref/iscontinuous.html
% clockChannel=s.addClockConnection('external',[deviceID,'/PFI12'],'ScanClock');
% %add Incoming clock connection https://de.mathworks.com/help/daq/clock-connections.html
lh = addlistener(s,'ErrorOccurred',@(src,event) disp(getReport(event.Error))); %will display otherwise invisible NIDAQ errors
s.queueOutputData(signal_min(1:Fs)); %put in the first data chunk, this is how many data points are in one second, https://de.mathworks.com/help/daq/ref/queueoutputdata.html % here put in our signal? 
all_fills=length(signal_min)/Fs-1; %Calc how many 1s-chunks your data contains. Total lenght of signal in data points divided by the sampling rates equals the number of 1 second data chunks Note: -1 because the first second already got queued
%all_fills=sig_dur-1; %Calc how many 1s-chunks your data contains. Note: -1 because the first second already got queued
n_bufferfill=0; %counter for how many chunks you already sent, will be updated by 1 each time the while loop is done 
disp('starting stimulation...')
secondCount=length(signal_min)/Fs; % shows how long in seconds the signal is, will be -1 each time the while loop is done, shows how many 1s chunks are left 
%secondCount=sig_dur; %to show time left
s.startBackground();        %starts the session, the NiDAQ will now generate data while depleting its buffer
running=true;
while running
    pause(0.05) %Essential pause. Necessary to update session in background. s.ScansQueued will not update otherwise
    if n_bufferfill<all_fills %as long as chunks of your data need sending ...
        if s.ScansQueued<Fs/2  %When only 500ms worth of data are left in the buffer,...
            n_bufferfill=n_bufferfill+1; %update your chunk-counter
            s.queueOutputData(signal_min(n_bufferfill*Fs+1:(n_bufferfill+1)*Fs)); %... put in the next chunk, calculates starting and end point of the next chunk 
            secondCount = secondCount - 1; %update time left in seconds 
            mle=floor(secondCount/60); %Conversion to minutes left, floor rounds to nearest integer 
            sle=mod(secondCount,60); %Conversion to seconds left, mod shows the remainder, e.g.: 123/60 = 2*60 and 3 seconds left while mle would show the 2 minutes 
            disp(sprintf('%d minute(s) and %d second(s) left', mle, sle));
        end
    elseif n_bufferfill==all_fills %If all chunks of your date were send 
        s.queueOutputData(zeros(Fs,1)); %fill zeros into buffer to avoid it running empty
        n_bufferfill=all_fills+1; %increase chunk counter so both if-clauses can't occur again
    elseif n_bufferfill>all_fills && s.ScansQueued<1e4 %If all chunks were sent, the extra zeros were put in and the buffer contains 1s 
        running=~running; %Session is done, not equal symbol  
    end
end
s.stop(); %kill session to avoid errors. stops s.startBackgound()


