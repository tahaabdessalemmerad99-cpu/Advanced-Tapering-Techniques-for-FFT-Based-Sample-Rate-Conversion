%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Python_like_Figure_AUDIO_ALL_Heatmap_WINDOWS_FINAL.m
% Code02
% FINAL VERSION
% - Real audio input
% - Fs_in extracted automatically from audio
% - Original paper tapers preserved
% - True Chebyshev implementation
% - True DDC implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUDIO INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fileName, pathName] = uigetfile('*.wav');

[file, Fs_in] = audioread( ...
    fullfile(pathName,fileName));

x = file(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs_out = 96000;

input_len = min(length(x),400000);

output_len = round(Fs_out / Fs_in * input_len);

x = x(1:input_len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = x ./ max(abs(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TAPER PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L =1212;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME AXIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_in = (0:input_len-1)' / Fs_in;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE TAPERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taper_naive = ones(input_len,1);

[taper_cosine,~] = ...
    get_taper('cosine',input_len,L);

[taper_hann,~] = ...
    get_taper('hann',input_len,L);

[taper_blackman,~] = ...
    get_taper('blackman',input_len,L);

[taper_cheb,~] = ...
    get_taper({'cheb',100},input_len,L);

[taper_ddc,~] = ...
    get_taper({'ddc',150},input_len,L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT SRC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_naive = fft_resample( ...
    x,taper_naive,output_len);

output_cosine = fft_resample( ...
    x,taper_cosine,output_len);

output_hann = fft_resample( ...
    x,taper_hann,output_len);

output_blackman = fft_resample( ...
    x,taper_blackman,output_len);

output_cheb = fft_resample( ...
    x,taper_cheb,output_len);

output_ddc = fft_resample( ...
    x,taper_ddc,output_len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY REGION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start = 0;
stop  = 5;
pad   = 0.05;

mask_in = ...
    (time_in >= start-pad) & ...
    (time_in < stop+pad);

time_out = (0:output_len-1)' / Fs_out;

mask_out = ...
    (time_out >= start-pad) & ...
    (time_out < stop+pad);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STFT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfft_in = 4096;
hop_in  = nfft_in/4;
win_in  = chebwin(nfft_in,150);

nfft_out = 4096;
hop_out  = nfft_out/4;
win_out  = chebwin(nfft_out,150);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w');

set(gcf,'Position',[100 100 1800 850]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) INPUT
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

hi = max(chunks,[],2);

t_down = ...
    (0:nBlocks-1)' * decimate / Fs_in;

fill([t_down; flipud(t_down)], ...
     [low; flipud(hi)], ...
     [0.2 0.2 0.2], ...
     'EdgeColor','none');

xlim([0 5]);

ylim([-1 1]);

title('(a)');

xlabel('Time (s)');

ylabel('$x(n)$', ...
    'Interpreter','latex');

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

xlim([start stop]);

ylim([20 22.5]);

caxis([-150 0]);

colormap(gca,parula);

hold on;

patch([start stop stop start], ...
      [Fs_out/2000 Fs_out/2000 ...
       Fs_in/2000 Fs_in/2000], ...
      [0.85 0.85 0.85], ...
      'EdgeColor','none');

plot([start stop], ...
     [Fs_out/2000 Fs_in/2000], ...
     'r-', ...
     'LineWidth',0.7);

plot([stop start], ...
     [Fs_out/2000 Fs_in/2000], ...
     'r-', ...
     'LineWidth',0.7);

uistack(findobj(gca,'Type','image'),'top');

title('(b)');

xlabel('Time (s)');

ylabel('Frequency (kHz)');

set(gca,'FontSize',11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
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
% OUTPUT SPECTROGRAMS
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

    ylim([20 22.5]);

    caxis([-150 0]);

    colormap(gca,parula);

    hold on;

    yline(Fs_in/2000, ...
        'w-', ...
        'LineWidth',0.8);

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

cb.TickLabels = ...
    {'-150','-100','-50','0 dB'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL TITLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sgtitle( ...
    ' SRC comparison (Audio)', ...
    'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = fft_resample(input,taper,output_len)

input = input(:);

input_len = length(input);

ratio = output_len / input_len;

input_f = fft(input);

if ratio >= 1

    input_f = input_f .* taper;

    output_f = zeros(output_len,1);

    pos_bins = floor((input_len+1)/2);

    neg_bins = floor((input_len-1)/2);

    output_f(1:pos_bins) = ...
        input_f(1:pos_bins);

    output_f(end-neg_bins+1:end) = ...
        input_f(end-neg_bins+1:end);

else

    output_f = zeros(output_len,1);

    pos_bins = floor((output_len+1)/2);

    neg_bins = floor((output_len-1)/2);

    output_f(1:pos_bins) = ...
        input_f(1:pos_bins);

    output_f(end-neg_bins+1:end) = ...
        input_f(end-neg_bins+1:end);

    output_f = output_f .* taper(1:output_len);

end

y = ratio * real(ifft(output_f));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [taper,meta] = ...
    get_taper(taper_spec,M,L)

taper = zeros(M,1);

meta = struct();

end_pos = floor((M+1)/2);

start_pos = end_pos - L + 1;

taper(1:end_pos) = 1;

if ~(ischar(taper_spec) && ...
        strcmpi(taper_spec,'box'))

    [transition,meta] = ...
        get_taper_transition( ...
        taper_spec,M,L);

    taper(start_pos:end_pos) = ...
        transition;

end

taper(floor((M+2)/2)+1:end) = ...
    taper(floor((M+1)/2):-1:2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [transition,meta] = ...
    get_taper_transition(win_spec,M,L)

[S,meta] = ...
    get_window(win_spec,M,L);

transition = flipud(cumsum(S(:)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S,meta] = ...
    get_window(win_spec,M,L)

meta = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COSINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(win_spec)

    switch lower(win_spec)

        case 'cosine'

            n = (0:L-1)';

            S = sin(pi*n/(L-1));

        case 'hann'

            n = (0:L-1)';

            S = 0.5 - ...
                0.5*cos(2*pi*n/(L-1));

        case 'blackman'

            n = (0:L-1)';

            S = ...
                0.42 ...
              - 0.5*cos(2*pi*n/(L-1)) ...
              + 0.08*cos(4*pi*n/(L-1));

        otherwise

            error('Unknown window.');

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEB / DDC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif iscell(win_spec)

    type = lower(win_spec{1});

    R = win_spec{2};

    switch type

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TRUE CHEBYSHEV
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        case 'cheb'

            order_m = L-1;

            sidelobes_R = 10^(R/20);

            x0 = ...
                cosh(acosh(sidelobes_R)/order_m);

            k = (0:order_m)';

            psi = x0 * cos(pi*k/L);

            T = zeros(size(psi));

            idx1 = abs(psi)<=1;

            idx2 = abs(psi)>1;

            T(idx1) = ...
                cos(order_m * acos(psi(idx1)));

            T(idx2) = ...
                cosh(order_m * acosh(abs(psi(idx2))));

            W = real(ifft(T));

            W = fftshift(W);

            W = abs(W);

            S = W(:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TRUE DDC
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        case 'ddc'

            alpha = ...
                binary_search_alpha(R,L);

            n = (0:L-1)'/(L-1);

            S = ...
                exp(-alpha*n.^2) ...
                .* ...
                exp(-alpha*(1-n).^2);

        otherwise

            error('Unknown tuple window.');

    end

else

    error('Invalid window specification.');

end

S = S(:);

S = S ./ sum(S);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function alpha = ...
    binary_search_alpha(R,L)

target = 10^(-R/20);

low = 0;

high = 200;

for k = 1:80

    alpha = (low+high)/2;

    n = (0:L-1)'/(L-1);

    w = ...
        exp(-alpha*n.^2) ...
        .* ...
        exp(-alpha*(1-n).^2);

    val = min(w);

    if val > target

        low = alpha;

    else

        high = alpha;

    end

end

alpha = (low+high)/2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S,F,T] = ...
    stft_nocenter(x,fs,nfft,hop,win)

x = x(:);

win = win(:);

n_frames = ...
    floor((length(x)-nfft)/hop)+1;

n_freq = floor(nfft/2)+1;

S = zeros(n_freq,n_frames);

for k = 1:n_frames

    idx = ...
        (1:nfft) + (k-1)*hop;

    frame = x(idx).*win;

    X = fft(frame,nfft);

    S(:,k) = X(1:n_freq);

end

F = (0:n_freq-1)' * (fs/nfft);

T = ((0:n_frames-1)*hop)/fs;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = amp_to_db(S)

A = 20*log10(abs(S)+1e-12);

A(A < -150) = -150;

end
