clc;
clearvars;
close all;

%% fit_levels_LG.m
% FINAL helper script for reproducing the department's historical
% Windows-volume calibration fit.
%
% This script is NOT required for routine participant calibration because
% adjust_volume_LG_A_final.m already contains the fitted coefficients.
%
% Run this script only if you want to inspect or reproduce the historical
% Windows-volume transfer curve.
%
% IMPORTANT:
% These measurements apply only to the audio hardware/output configuration
% in which they were acquired. If the computer, audio interface, headphones,
% driver processing, or other parts of the playback chain have changed, the
% physical output should be measured again.

% Historical acoustic level assigned to Windows slider 100.
dBA_100 = 70.5;

% Historical measurements:
% [Windows slider position, measured dB(A)]
recorded_levels = [ ...
      1   45.1
      2   54.5
      3   60.2
      4   64.3
      5   67.6
      7   72.6
     10   78.0
     15   84.2
     20   88.5
     25   91.0
     30   94.8
     40   99.2
     50  102.7
     60  105.5
     70  107.9
     80  110.0
     90  111.8
    100  113.4
];

% Apply the historical correction so slider 100 corresponds to dBA_100.
correction = recorded_levels(end,2) - dBA_100;

x = recorded_levels(:,1);
y = recorded_levels(:,2) - correction;

% Fit:
%       level = a + b*ln(Windows slider)
%
% Ordinary least squares is used, so no Curve Fitting Toolbox is required.
X = [ones(size(x)) log(x)];
beta = X \ y;

a = beta(1);
b = beta(2);

fprintf('Fitted Windows-volume relationship:\n');
fprintf('  level = %.4f + %.4f * ln(slider)\n\n', a, b);

fprintf('Coefficients used in adjust_volume_LG_A_final.m:\n');
fprintf('  level = 0.9236 + 15.01 * ln(slider)\n\n');

fprintf('Historical reference at Windows 100: %.1f dB(A)\n', dBA_100);
fprintf('Applied correction to raw measurements: %.1f dB\n\n', correction);

% Plot corrected historical measurements and fitted curve.
x1 = linspace(1,100,500)';
y1 = a + b*log(x1);

figure;
hold on;
plot(x,y,'ok','MarkerFaceColor','none');
plot(x1,y1,'LineWidth',1.5);
legend({'Measurements','Logarithmic fit'},'Location','best');
xlabel('Windows volume slider position');
ylabel('Corrected acoustic level [dB(A)]');
title('Windows volume slider calibration');
grid on;
