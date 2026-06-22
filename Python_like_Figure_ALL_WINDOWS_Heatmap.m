%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Python_like_Figure_ALL_WINDOWS_Heatmap.m
%
% Faithful MATLAB translation of:
% example.py / taper.py / resamp.py
%
% WINDOWS:
% - No taper
% - Cosine
% - Hann
% - Blackman
% - Chebyshev
% - DDC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

rng(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs_in  = 8000;
Fs_out = 22050;
%Fs_out = 4000;
input_len  = 400000;
output_len = round(Fs_out / Fs_in * input_len);

L = round(0.10 * input_len / 2);

R = 150;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PYTHON-LIKE INPUT SIGNAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = 2*rand(input_len,1) - 1;

time_in = (0:input_len-1)' / Fs_in;

envelope = exp(-40 * mod(time_in,1));

x = x .* envelope;

x = x ./ max(abs(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESAMPLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NO TAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_naive = fft_resample( ...
    x, ...
    ones(input_len,1), ...
    output_len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COSINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taper_cosine = get_taper('cosine', input_len, L);

output_cosine = fft_resample( ...
    x, ...
    taper_cosine, ...
    output_len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taper_hann = get_taper('hann', input_len, L);

output_hann = fft_resample( ...
    x, ...
    taper_hann, ...
    output_len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLACKMAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taper_blackman = get_taper('blackman', input_len, L);

output_blackman = fft_resample( ...
    x, ...
    taper_blackman, ...
    output_len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEBYSHEV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taper_cheb = get_taper({'chebwin',R}, input_len, L);

output_cheb = fft_resample( ...
    x, ...
    taper_cheb, ...
    output_len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DDC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taper_ddc = get_taper({'ddc',R}, input_len, L);

output_ddc = fft_resample( ...
    x, ...
    taper_ddc, ...
    output_len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING WINDOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start = 1.0;

stop  = 3.0;

pad   = 0.05;

mask_in = ...
    (time_in >= start-pad) & ...
    (time_in < stop+pad);

time_out = (0:output_len-1)' / Fs_out;

mask_out = ...
    (time_out >= start-pad) & ...
    (time_out < stop+pad);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STFT SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfft_in  = 512;
hop_in   = nfft_in / 4;
win_in   = chebwin(nfft_in,160);

nfft_out = 1024;
hop_out  = nfft_out / 4;
win_out  = chebwin(nfft_out,160);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w');

set(gcf,'Position',[100 100 1800 850]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) INPUT WAVEFORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,4,1);

input_half = x(1:min(length(x),5*Fs_in));

decimate = 20;

nBlocks = floor(length(input_half)/decimate);

chunks = reshape( ...
    input_half(1:nBlocks*decimate), ...
    decimate, ...
    nBlocks).';

low = min(chunks,[],2);

hi  = max(chunks,[],2);

t_down = (0:nBlocks-1)' * decimate / Fs_in;

fill([t_down; flipud(t_down)], ...
     [low; flipud(hi)], ...
     [0.2 0.2 0.2], ...
     'EdgeColor','none');

xlim([0 5]);

ylim([-1 1]);

title('(a)');

xlabel('Time (s)');

ylabel('$x(n)$','Interpreter','latex');

set(gca,'FontSize',11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (b) INPUT SPECTROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,4,5);

x_in_seg = x(mask_in);

[S1,F1,T1] = stft_nocenter( ...
    x_in_seg, ...
    Fs_in, ...
    nfft_in, ...
    hop_in, ...
    win_in);

T1 = T1 + start;

SdB1 = amp_to_db(S1);

SdB1 = SdB1 - max(SdB1(:));

SdB1(SdB1 < -150) = -150;

imagesc(T1,F1/1000,SdB1);

axis xy;

xlim([start-0.02 stop+0.02]);

ylim([0 Fs_out/2000]);

caxis([-150 0]);

colormap(gca,parula);

hold on;

patch([start stop stop start], ...
      [Fs_in/2000 Fs_in/2000 Fs_out/2000 Fs_out/2000], ...
      [0.85 0.85 0.85], ...
      'EdgeColor','none', ...
      'FaceAlpha',1.0);

plot([start stop], ...
     [Fs_in/2000 Fs_out/2000], ...
     'r-', ...
     'LineWidth',0.7);

plot([stop start], ...
     [Fs_in/2000 Fs_out/2000], ...
     'r-', ...
     'LineWidth',0.7);

uistack(findobj(gca,'Type','image'),'top');

title('(b)');

xlabel('Time (s)');

ylabel('Frequency (kHz)');

set(gca,'FontSize',11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT SIGNALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

signals = {
    output_naive      '(c) No taper';
    output_cosine     '(d) Cosine';
    output_hann       '(e) Hann';
    output_blackman   '(f) Blackman';
    output_cheb       '(g) Chebyshev';
    output_ddc        '(h) DDC';
};

positions = [2 3 4 6 7 8];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:6

    subplot(2,4,positions(k));

    y = signals{k,1};

    ttl = signals{k,2};

    y_seg = y(mask_out);

    [S,F,T] = stft_nocenter( ...
        y_seg, ...
        Fs_out, ...
        nfft_out, ...
        hop_out, ...
        win_out);

    T = T + start;

    SdB = amp_to_db(S);

    SdB = SdB - max(SdB(:));

    SdB(SdB < -150) = -150;

    imagesc(T,F/1000,SdB);

    axis xy;

    xlim([start stop]);

    ylim([0 Fs_out/2000]);

    caxis([-150 0]);

    colormap(gca,parula);

    hold on;

    yline(Fs_in/1000, ...
        'w-', ...
        'LineWidth',0.7);

    plot([1.15 1.35], ...
         [6.5 4.1], ...
         'w-', ...
         'LineWidth',2);

    plot([1.35 1.30], ...
         [4.1 4.25], ...
         'w-', ...
         'LineWidth',2);

    plot([1.35 1.28], ...
         [4.1 4.15], ...
         'w-', ...
         'LineWidth',2);

    title(ttl);

    xlabel('Time (s)');

    ylabel('Frequency (kHz)');

    set(gca,'FontSize',11);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLORBAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cb = colorbar;

cb.Position = [0.93 0.12 0.02 0.75];

cb.Ticks = [-150 -100 -50 0];

cb.TickLabels = {'-150','-100','-50','0 dB'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN TITLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sgtitle( ...
    ' SRC comparison', ...
    'FontWeight', ...
    'bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = fft_resample(input, taper, output_len)

input = input(:);

ratio = output_len / length(input);

input_f = fft(input);

if ratio >= 1

    if ~isempty(taper)
        input_f = input_f .* taper;
    end

    output_f = zeros(output_len,1,'like',input_f);

    pos_bins = floor((length(input)+1)/2);

    neg_bins = floor((length(input)-1)/2);

    output_f(1:pos_bins) = input_f(1:pos_bins);

    output_f(end-neg_bins+1:end) = ...
        input_f(end-neg_bins+1:end);

else

    output_f = zeros(output_len,1,'like',input_f);

    pos_bins = floor((output_len+1)/2);

    neg_bins = floor((output_len-1)/2);

    output_f(1:pos_bins) = input_f(1:pos_bins);

    output_f(end-neg_bins+1:end) = ...
        input_f(end-neg_bins+1:end);

    if ~isempty(taper)
        output_f = output_f .* taper;
    end

end

y = ratio * real(ifft(output_f));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W, meta] = get_taper(taper_spec, M, L)

W = zeros(M,1);

meta = struct();

end_pos = floor((M+1)/2);

start_pos = end_pos - L + 1;

W(1:end_pos) = 1;

if ~(ischar(taper_spec) && strcmpi(taper_spec,'box'))

    [transition, meta] = ...
        get_taper_transition(taper_spec, M, L);

    W(start_pos:end_pos) = transition;

end

W(floor((M+2)/2)+1:end) = ...
    W(floor((M+1)/2):-1:2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [transition, meta] = ...
    get_taper_transition(win_spec, M, L)

[S, meta] = get_window(win_spec, M, L);

meta.window = S;

transition = flipud(cumsum(S(:)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S, meta] = get_window(win_spec, M, L)

meta = struct();

if ischar(win_spec) || isstring(win_spec)

    switch lower(char(win_spec))

        case 'cosine'

            n = (0:L-1)';

            S = sin(pi*n/(L-1));

        case 'hann'

            S = hann(L,'symmetric');

        case 'blackman'

            S = blackman(L,'symmetric');

        otherwise

            error('Unknown window');

    end

elseif iscell(win_spec)

    type = lower(win_spec{1});

    switch type

        case 'chebwin'

            sidelobes_dB = win_spec{2};

            sidelobes_R = db2amp(sidelobes_dB);

            order_m = L - 1;

            x0 = cosh(acosh(sidelobes_R)/order_m);

            Lr = floor(L/2)+1;

            t = (0:Lr-1)' * pi / L;

            x = x0 * cos(t);

            s = complex(zeros(Lr,1));

            mask = abs(x)<=1;

            s(mask) = cos(order_m * acos(x(mask)));

            mask = x>1;

            s(mask) = cosh(order_m * acosh(x(mask)));

            mask = x<-1;

            s(mask) = (-1)^order_m * ...
                cosh(order_m * acosh(-x(mask)));

            if mod(L,2)==0
                s = s .* exp(1j*t);
            end

            if mod(L,2)==0
                fullspec = [s; conj(s(end-1:-1:2))];
            else
                fullspec = [s; conj(s(end:-1:2))];
            end

            S = fftshift(real(ifft(fullspec,L)));

            S = S(:);

            S = S / sum(S);

        case 'ddc'

            sidelobes_dB = win_spec{2};

            sidelobes_R = db2amp(sidelobes_dB);

            order_m = L - 1;

            x0 = binary_search( ...
                @(z) ddc_DQ(z,order_m,M,L), ...
                [1 2], ...
                sidelobes_R);

            alpha = ...
                (1 + sqrt(x0^2 - 1)/x0)/2;

            Lr = floor(L/2)+1;

            t = (0:Lr-1)' * pi / L;

            x = x0 * cos(t);

            s = complex(zeros(Lr,1));

            mask = abs(x)<=1;

            s(mask) = ...
                alpha * cos(order_m * acos(x(mask))) ...
              - (1-alpha) * cos((order_m-2) * acos(x(mask)));

            mask = x>1;

            s(mask) = ...
                alpha * cosh(order_m * acosh(x(mask))) ...
              - (1-alpha) * cosh((order_m-2) * acosh(x(mask)));

            mask = x<-1;

            s(mask) = (-1)^order_m * ...
                (alpha * cosh(order_m * acosh(-x(mask))) ...
              - (1-alpha) * cosh((order_m-2) * acosh(-x(mask))));

            if mod(L,2)==0
                s = s .* exp(1j*t);
            end

            if mod(L,2)==0
                fullspec = [s; conj(s(end-1:-1:2))];
            else
                fullspec = [s; conj(s(end:-1:2))];
            end

            S = fftshift(real(ifft(fullspec,L)));

            S = S(:);

            S = S / sum(S);

        otherwise

            error('Unknown tuple window.');

    end

else

    error('Invalid window specification.');

end

S = S(:);

S = S / sum(S);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = ddc_DQ(x0,order_m,M,L)

alpha = ...
    (1 + sqrt(x0.^2 - 1)./x0)/2;

Q = ...
    alpha * cosh(order_m * acosh(x0)) ...
  - (1-alpha) * cosh((order_m-2) * acosh(x0));

val = (M-L) * Q;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = binary_search(f,x_lim,y)

delta = (x_lim(2)-x_lim(1))/2;

x = x_lim(1);

while delta > 1e-20

    if f(x + delta) <= y
        x = x + delta;
    end

    delta = delta / 2;

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = db2amp(d)

y = 10.^(d/20);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S,F,T] = stft_nocenter(x,fs,nfft,hop,win)

x = x(:);

win = win(:);

n_frames = floor((length(x)-nfft)/hop)+1;

n_freq = floor(nfft/2)+1;

S = zeros(n_freq,n_frames);

for k = 1:n_frames

    idx = (1:nfft) + (k-1)*hop;

    frame = x(idx) .* win;

    X = fft(frame,nfft);

    S(:,k) = X(1:n_freq);

end

F = (0:n_freq-1)' * (fs/nfft);

T = ((0:n_frames-1)*hop)/fs;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = amp_to_db(S)

A = 20*log10(abs(S)+1e-10);

A(A < -150) = -150;

end