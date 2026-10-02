function [EEG_post, figPSD] = clean_line_noise(EEG_pre, varargin)
% LG_CLEAN_LINE_NOISE  Remove power-line noise (CleanLine or notch) and plot pre/post PSD.
%
%   [EEG_post, figPSD] = lg_clean_line_noise(EEG_pre, 'Name',Value, ...)
%
% Inputs
%   EEG_pre        : EEGLAB EEG struct (continuous)
%
% Name-Value pairs (all optional)
%   'linefreq'     : fundamental line frequency (default 50)
%   'useCleanline' : try pop_cleanline if available (default true)
%   'bandwidthHz'  : CleanLine half-bandwidth in Hz (default 2)
%   'headerVhdr'   : path to the .vhdr to auto-skip if an online notch is detected (default '')
%   'chanlist'     : channel indices to clean (default 1:EEG.nbchan)
%   'maxHarmonics' : up to N harmonics (CleanLine target list) (default 4)
%   'outPng'       : filename to save the pre/post PSD report (default '')
%   'plotFmax'     : max frequency shown in PSD plot (default min(100, fs/2))
%
% Outputs
%   EEG_post : EEG after line-noise cleanup
%   figPSD   : handle to PSD comparison figure (or [] if plotting is skipped)
%
% Requirements
%   CleanLine plugin (optional). If not present, the function falls back to a
%   narrow band-stop notch at ±1 Hz around the fundamental (and the 2nd harmonic).
%
% Miguel-friendly, 2025 :)

% ---------- defaults ----------
p = inputParser;
p.addParameter('linefreq',     50,     @(x)isnumeric(x)&&isscalar(x));
p.addParameter('useCleanline', true,   @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('bandwidthHz',  2,      @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('headerVhdr',   '',     @(x)ischar(x)||isstring(x));
p.addParameter('chanlist',     [],     @(x)isnumeric(x)||islogical(x)||isempty(x));
p.addParameter('maxHarmonics', 4,      @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('outPng',       '',     @(x)ischar(x)||isstring(x));
p.addParameter('plotFmax',     [],     @(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;

EEG_post  = EEG_pre;   % copy
figPSD    = [];

if isempty(opt.chanlist)
    opt.chanlist = 1:EEG_pre.nbchan;
end
fs = EEG_pre.srate;

% --- skip if header indicates an online notch was used ---
if ~isempty(opt.headerVhdr) && header_has_notch(opt.headerVhdr)
    fprintf('LG_CLEAN_LINE_NOISE: Online notch detected in header -> skipping offline cleanup.\n');
else
    fprintf('LG_CLEAN_LINE_NOISE: Removing line noise at %d Hz...\n', opt.linefreq);

    % --- try CleanLine first ---
    cleaned = false;
    if opt.useCleanline && exist('pop_cleanline','file') == 2
        maxH = min(fs/2, opt.maxHarmonics*opt.linefreq);
        tgt  = opt.linefreq:opt.linefreq:maxH;
        try
            EEG_post = pop_cleanline(EEG_post, ...
                'bandwidth',     opt.bandwidthHz, ...
                'chanlist',      opt.chanlist, ...
                'computepower',  1, ...
                'linefreqs',     tgt, ...
                'normSpectrum',  0, ...
                'p',             0.01, ...
                'pad',           2, ...
                'plotfigures',   0, ...
                'scanforlines',  1, ...
                'sigtype',       'Channels', ...
                'tau',           100, ...
                'verb',          0);
            cleaned = true;
        catch
            warning('CleanLine failed. Falling back to notch.');
        end
    end

    % --- fallback: narrow notch at fundamental ±1 Hz (and 2nd harmonic) ---
    if ~cleaned
        span1 = [opt.linefreq-1, opt.linefreq+1];
        EEG_post = pop_eegfiltnew(EEG_post, 'locutoff', span1(1), 'hicutoff', span1(2), ...
                                  'revfilt', 1, 'plotfreqz', 0);
        if fs/2 >= 2*opt.linefreq   % optional 2nd harmonic
            span2 = [2*opt.linefreq-1, 2*opt.linefreq+1];
            EEG_post = pop_eegfiltnew(EEG_post, 'locutoff', span2(1), 'hicutoff', span2(2), ...
                                      'revfilt', 1, 'plotfreqz', 0);
        end
    end
    EEG_post = eeg_checkset(EEG_post);
end

% --- PSD comparison plot (pre vs post) ---
try
    figPSD = plot_line_noise_psd(EEG_pre, EEG_post, opt.linefreq, opt.outPng, opt.plotFmax);
catch
    warning('PSD plotting failed.');
end

end % main function


% ===================== helpers =====================

function hadNotch = header_has_notch(vhdrFullPath)
% Return true if BrainVision header suggests a notch was applied online.
hadNotch = false;
try
    txt = fileread(vhdrFullPath);
    hadNotch = ~isempty(regexpi(txt,'Notch.*(On|Yes|True)','once')) || ...
               ~isempty(regexpi(txt,'Filter.*Notch.*(On|Yes|True)','once'));
catch
    hadNotch = false;
end
end


function fig = plot_line_noise_psd(EEG_pre, EEG_post, linefreq, out_png, plotFmax)
% Overlay median PSD (pre vs post), scalp map of per-channel reduction @ fundamental,
% and sorted per-channel bar chart.

assert(EEG_pre.nbchan == EEG_post.nbchan, 'Channel mismatch pre/post.');
assert(EEG_pre.srate  == EEG_post.srate,  'Sampling rate mismatch pre/post.');

fs   = EEG_pre.srate;
if isempty(plotFmax), plotFmax = min(100, fs/2); end
frng = [1 plotFmax];

% Spectra (dB) per channel
[spec_pre, freqs]  = spectopo(double(EEG_pre.data),  0, fs, 'plot','off', 'freqrange', frng);
[spec_post, ~]     = spectopo(double(EEG_post.data), 0, fs, 'plot','off', 'freqrange', frng);

med_pre  = median(spec_pre, 1);
med_post = median(spec_post, 1);

% indices of fundamental + harmonics
harm = linefreq:linefreq:plotFmax;
idxh = arrayfun(@(f) find(abs(freqs - f)==min(abs(freqs - f)),1,'first'), harm);

% dB reduction at fundamental
idx0          = idxh(1);
chan_drop_db  = spec_pre(:,idx0) - spec_post(:,idx0);   % + = reduction
med_drop_db   = median(chan_drop_db);

% ---- figure ----
fig = figure('Units','normalized','Position',[0 0 1 1],'Color','w');
t   = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

% 1) PSD overlay
ax1 = nexttile(t,[1 2]); hold(ax1,'on');
plot(ax1, freqs, med_pre,  'LineWidth',1.6);
plot(ax1, freqs, med_post, 'LineWidth',1.6);
yline(ax1,0,'k-'); grid(ax1,'on');
xlabel(ax1,'Frequency (Hz)'); ylabel(ax1,'Power (dB)');
title(ax1, sprintf('Median PSD across channels (pre vs post)  |  Δ @ %d Hz = %.1f dB', linefreq, med_drop_db));
legend(ax1, {'Pre-clean','Post-clean'}, 'Location','northeast','Box','off');
yl = ylim(ax1);
for k=1:numel(harm)
    xline(ax1, harm(k), ':', 'Color',[0.4 0 0], 'LineWidth',0.8, 'HandleVisibility','off');
end
ylim(ax1, yl);

% 2) scalp map (reduction @ fundamental)
ax2 = nexttile(t,3); axes(ax2);
ws = warning('off','MATLAB:ui:Axes:UnsupportedConveniencePositionSet');
topoplot(chan_drop_db, EEG_post.chanlocs, ...
    'electrodes','on', 'emarker',{'.','k',6,1}, ...
    'numcontour',6, 'gridscale',120);
warning(ws);
cb = colorbar(ax2,'eastoutside');
cb.Label.String = sprintf('dB reduction @ %d Hz (pre - post)', linefreq);
title(ax2, 'Reduction map (higher is better)');

% 3) per-channel sorted bar
ax3 = nexttile(t,4);
[sorted_drop,~] = sort(chan_drop_db,'descend');
bar(ax3, sorted_drop, 'FaceColor',[0.2 0.5 0.9]);
yline(ax3,0,'k-'); grid(ax3,'on');
xlabel(ax3,'Channels (sorted)'); ylabel(ax3,'dB reduction');
title(ax3,'Per-channel reduction @ fundamental');
xlim(ax3,[0 numel(sorted_drop)+1]);

sgtitle(t, sprintf('Line-noise check (Pre vs Post)  |  %d Hz (and harmonics)', linefreq), ...
    'FontWeight','bold');

% save if requested
if nargin>=5 && ~isempty(out_png)
    try
        exportgraphics(fig, out_png, 'Resolution', 200);
    catch
        saveas(fig, out_png);
    end
end
end
