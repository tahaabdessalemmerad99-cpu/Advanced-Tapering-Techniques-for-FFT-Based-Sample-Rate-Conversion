%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare_tapers_impulse_response.m
%
% Impulse response comparison for:
% - No taper
% - Cosine
% - Hann
% - Blackman
% - Chebyshev
% - DDC
%
% Uses EXACT SAME FUNCTIONS as final FFT SRC code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 44100;

M = 192000;

L = 1102;

R = 150;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPULSE INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = zeros(N,1);

x(N/2) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TAPER DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tapers = {
    'No taper',  'box';
    'Cosine',    'cosine';
    'Hann',      'hann';
    'Blackman',  'blackman';
    'Chebyshev', {'chebwin',95.8};
    'DDC',       {'ddc',150};
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cols = lines(size(tapers,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME AXIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = ((0:M-1) - M/2)/M * 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w');

set(gcf,'Position',[100 100 1400 800]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP OVER TAPERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:size(tapers,1)

    name = tapers{k,1};

    spec = tapers{k,2};

    fprintf('\nProcessing: %s\n', name);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CREATE TAPER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ischar(spec) && strcmpi(spec,'box')

        taper = ones(N,1);

    else

        taper = get_taper(spec, N, L);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FFT RESAMPLE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    y = fft_resample(x, taper, M);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NORMALIZATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    y = y ./ max(abs(y));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SUBPLOT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subplot(2,3,k);

    plot(time, amp2db_custom(y), ...
        'LineWidth', 1.2, ...
        'Color', cols(k,:));

    grid on;

    xlim([-15 15]);

    ylim([-220 10]);

    xlabel('Time (ms)');

    ylabel('Magnitude (dB)');

    title(name);

    set(gca,'FontSize',11);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL TITLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sgtitle('Impulse Response Comparison of Tapers');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT RESAMPLE
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
% GET TAPER
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
% GET TRANSITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [transition, meta] = ...
    get_taper_transition(win_spec, M, L)

[S, meta] = get_window(win_spec, M, L);

meta.window = S;

transition = flipud(cumsum(S(:)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET WINDOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S, meta] = get_window(win_spec, M, L)

meta = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STANDARD WINDOWS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIAL WINDOWS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif iscell(win_spec)

    type = lower(win_spec{1});

    switch type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEBYSHEV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DDC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% DDC DQ
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
% BINARY SEARCH
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
% DB TO AMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = db2amp(d)

y = 10.^(d/20);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AMP TO DB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = amp2db_custom(x)

y = 20*log10(abs(x)+1e-15);

end