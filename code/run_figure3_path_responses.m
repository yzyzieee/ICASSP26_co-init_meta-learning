% Copyright (c) 2026 Ziyi Yang (ziyi016@e.ntu.edu.sg).
% Released for the ICASSP 2026 MAML co-initialization experiments.
%
% Reproduce Fig. 3: measured primary and secondary path responses.
%
% The highlighted rows correspond to the two three-task training sets used
% in the paper. All available PANDAR responses are shown in gray.

clear; clc; close all;

code_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(code_dir);
data_dir = fullfile(root_dir, 'data');
out_dir  = fullfile(root_dir, 'results');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

pairs_file = fullfile(data_dir, 'ANC_PriSec_pairs.mat');

nfft = 4096;
f_band = [100 2000];
f_plot_max = 2000;
eps0 = 1e-12;

% Training set A: diverse path set. Training set B: compact path set.
rowsA = [32 37 43];
rowsB = [20 26 31];

S = load(pairs_file, 'out_pairs');
H_primary = S.out_pairs.H_primary;
H_secondary = S.out_pairs.H_secondary;
meta = S.out_pairs.meta_pairs;
Fs0 = S.out_pairs.Fs;

[magP, f] = response_bank(H_primary, meta.prim_col, Fs0, nfft, eps0);
[magS, ~] = response_bank(H_secondary, meta.sec_col, Fs0, nfft, eps0);

fig = figure('Color', 'w', 'Position', [80 80 1200 680]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile; hold(ax1, 'on'); grid(ax1, 'on');
plot_all_responses(ax1, f, magP, [0.72 0.72 0.72]);
plot_group(ax1, f, H_primary, meta.prim_col(rowsA), Fs0, nfft, eps0, [0.85 0.20 0.20], '-');
plot_group(ax1, f, H_primary, meta.prim_col(rowsB), Fs0, nfft, eps0, [0.20 0.45 0.85], '--');
title(ax1, '(a) Primary paths');
ylabel(ax1, 'Magnitude (dB)');
xlim(ax1, [0 f_plot_max]);
ylim(ax1, [-50 20]);
add_band(ax1, f_band);
add_path_legend(ax1, 'northeast');

ax2 = nexttile; hold(ax2, 'on'); grid(ax2, 'on');
plot_all_responses(ax2, f, magS, [0.72 0.72 0.72]);
plot_group(ax2, f, H_secondary, meta.sec_col(rowsA), Fs0, nfft, eps0, [0.85 0.20 0.20], '-');
plot_group(ax2, f, H_secondary, meta.sec_col(rowsB), Fs0, nfft, eps0, [0.20 0.45 0.85], '--');
title(ax2, '(b) Secondary paths');
xlabel(ax2, 'Frequency (Hz)');
ylabel(ax2, 'Magnitude (dB)');
xlim(ax2, [0 f_plot_max]);
ylim(ax2, [-50 20]);
add_band(ax2, f_band);
add_path_legend(ax2, 'southeast');

set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(fig, '-property', 'FontSize'), 'FontSize', 14);

save_figure(fig, fullfile(out_dir, 'Fig3_path_responses.png'));

function [mag_db, f] = response_bank(Hmat, meta_cols, Fs, nfft, eps0)
    valid = isfinite(meta_cols) & meta_cols >= 1 & meta_cols <= size(Hmat, 2);
    cols = unique(meta_cols(valid));
    mag_db = zeros(nfft, numel(cols));
    f = [];
    for k = 1:numel(cols)
        h = double(Hmat(:, cols(k)));
        [H, f] = freqz(h, 1, nfft, Fs);
        mag_db(:, k) = 20*log10(abs(H) + eps0);
    end
end

function plot_all_responses(ax, f, mag_db, color)
    for k = 1:size(mag_db, 2)
        plot(ax, f, mag_db(:, k), 'Color', color, 'LineWidth', 0.7);
    end
end

function plot_group(ax, f, Hmat, cols, Fs, nfft, eps0, color, line_style)
    cols = cols(:).';
    for k = 1:numel(cols)
        h = double(Hmat(:, cols(k)));
        H = freqz(h, 1, nfft, Fs);
        y = 20*log10(abs(H) + eps0);
        plot(ax, f, y, 'Color', color, 'LineWidth', 2.2, 'LineStyle', line_style);
    end
end

function add_band(ax, f_band)
    yl = ylim(ax);
    p = patch(ax, [f_band(1) f_band(2) f_band(2) f_band(1)], ...
        [yl(1) yl(1) yl(2) yl(2)], [0.9 0.95 1.0], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.12, 'HitTest', 'off');
    uistack(p, 'bottom');
end

function add_path_legend(ax, loc)
    h1 = plot(ax, NaN, NaN, 'Color', [0.72 0.72 0.72], 'LineWidth', 0.9);
    h2 = plot(ax, NaN, NaN, 'Color', [0.85 0.20 0.20], 'LineWidth', 2.4, 'LineStyle', '-');
    h3 = plot(ax, NaN, NaN, 'Color', [0.20 0.45 0.85], 'LineWidth', 2.4, 'LineStyle', '--');
    legend(ax, [h1 h2 h3], {'All responses', 'Training set A', 'Training set B'}, ...
        'Location', loc);
end

function save_figure(fig, filename)
    try
        exportgraphics(fig, filename, 'Resolution', 300);
    catch
        saveas(fig, filename);
    end
    fprintf('Saved %s\n', filename);
end
