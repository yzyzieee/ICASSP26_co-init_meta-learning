% Copyright (c) 2026 Ziyi Yang (ziyi016@e.ntu.edu.sg).
% Released for the ICASSP 2026 MAML co-initialization experiments.
%
% Reproduce Fig. 2: three-phase OSPM-FxLMS with path switches.
%
% This script uses the paired PANDAR primary/secondary paths and the
% pre-trained MAML co-initialization stored in ../data.

clear; clc; close all;

code_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(code_dir);
data_dir = fullfile(root_dir, 'data');
out_dir  = fullfile(root_dir, 'results');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

fs = 16000;
pairs_file = fullfile(data_dir, 'ANC_PriSec_pairs.mat');
init_file  = fullfile(data_dir, 'Wc_Sc_MAML_RWTH_Min_3.mat');

% These are the three evaluation rows used for the switch experiment.
task_rows = [10 25 38];

band_hz = [200 1000];
fir_ord = 512;

Lw = 512;
Ls = 512;
Lh = 512;
delay = 0;
snr_db = 30;
phase_sec = [60 60 60];
phase_len = phase_sec * fs;

alpha_power = 0.99;
c_aux = 1.0;
mu_w = 1e-2;
mu_s = 1e-3;
mu_h = 1e-3;
Tw = 10;
Ts = 10;
Th = 10;

init = load(init_file, 'Wc', 'Sc');
w_maml = pad_or_trim(init.Wc(:), Lw);
s_maml = pad_or_trim(init.Sc(:), Ls);

pairs = load(pairs_file, 'out_pairs');
Fs0 = pairs.out_pairs.Fs;
Hpri = pairs.out_pairs.H_primary;
Hsec = pairs.out_pairs.H_secondary;
meta = pairs.out_pairs.meta_pairs;

[P1, S1] = pick_pair(meta, Hpri, Hsec, task_rows(1), fs, Fs0);
[P2, S2] = pick_pair(meta, Hpri, Hsec, task_rows(2), fs, Fs0);
[P3, S3] = pick_pair(meta, Hpri, Hsec, task_rows(3), fs, Fs0);

p1 = pad_or_trim(P1(:), Lw);  s1 = pad_or_trim(S1(:), Ls);
p2 = pad_or_trim(P2(:), Lw);  s2 = pad_or_trim(S2(:), Ls);
p3 = pad_or_trim(P3(:), Lw);  s3 = pad_or_trim(S3(:), Ls);

bp = fir1(fir_ord, band_hz/(fs/2));
n_total = sum(phase_len);
ref = filter(bp, 1, randn(n_total, 1));
ref_noisy = add_awgn_measured(ref, snr_db);

idx1 = 1:phase_len(1);
idx2 = phase_len(1)+(1:phase_len(2));
idx3 = sum(phase_len(1:2))+(1:phase_len(3));

d1 = add_awgn_measured(filter(p1, 1, ref(idx1)), snr_db);
d2 = add_awgn_measured(filter(p2, 1, ref(idx2)), snr_db);
d3 = add_awgn_measured(filter(p3, 1, ref(idx3)), snr_db);

ref1 = ref_noisy(idx1);
ref2 = ref_noisy(idx2);
ref3 = ref_noisy(idx3);

% Match the final paper experiment where later phases have higher level.
gain2 = 10^(2.5/20);
gain3 = 10^(3.5/20);
ref2 = ref2 * gain2;  d2 = d2 * gain2;
ref3 = ref3 * gain3;  d3 = d3 * gain3;

rng(1234);
[eZ1, stZ1] = run_osmp_phase(ref1, d1, s1, ...
    zeros(Lw,1), zeros(Ls,1), zeros(Lh,1), ...
    mu_w, mu_s, mu_h, Tw, Ts, Th, alpha_power, c_aux, delay);

rng(1234);
[eM1, stM1] = run_osmp_phase(ref1, d1, s1, ...
    w_maml, s_maml, zeros(Lh,1), ...
    mu_w, mu_s, mu_h, Tw, Ts, Th, alpha_power, c_aux, delay);

rng(5678);
[eZ2, stZ2] = run_osmp_phase(ref2, d2, s2, ...
    stZ1.w, stZ1.s, stZ1.h, ...
    mu_w, mu_s, mu_h, Tw, Ts, Th, alpha_power, c_aux, delay);

rng(5678);
[eM2, stM2] = run_osmp_phase(ref2, d2, s2, ...
    w_maml, s_maml, zeros(Lh,1), ...
    mu_w, mu_s, mu_h, Tw, Ts, Th, alpha_power, c_aux, delay);

rng(9012);
[eZ3, stZ3] = run_osmp_phase(ref3, d3, s3, ...
    stZ2.w, stZ2.s, stZ2.h, ...
    mu_w, mu_s, mu_h, Tw, Ts, Th, alpha_power, c_aux, delay);

rng(9012);
[eM3, stM3] = run_osmp_phase(ref3, d3, s3, ...
    w_maml, s_maml, zeros(Lh,1), ...
    mu_w, mu_s, mu_h, Tw, Ts, Th, alpha_power, c_aux, delay);

e_original = [eZ1; eZ2; eZ3];
e_maml = [eM1; eM2; eM3];
d_all = [d1; d2; d3];
aux_original = [stZ1.aux; stZ2.aux; stZ3.aux];
aux_maml = [stM1.aux; stM2.aux; stM3.aux];

win = 4096;
ms_off = movmean(d_all.^2, win, 'Endpoints', 'shrink');
ms_original = movmean(e_original.^2, win, 'Endpoints', 'shrink');
ms_maml = movmean(e_maml.^2, win, 'Endpoints', 'shrink');
ms_aux_original = movmean(aux_original.^2, win, 'Endpoints', 'shrink');
ms_aux_maml = movmean(aux_maml.^2, win, 'Endpoints', 'shrink');

t = (0:numel(d_all)-1)/fs;
t_switch1 = phase_sec(1);
t_switch2 = phase_sec(1) + phase_sec(2);

fig = figure('Color', 'w', 'Position', [80 120 1500 360]);
hold on; grid on;
plot(t, 10*log10(ms_off + 1e-12), 'k', 'LineWidth', 1.6);
plot(t, 10*log10(ms_original + 1e-12), '--', 'LineWidth', 1.6, 'Color', [0.8500 0.3250 0.0980]);
plot(t, 10*log10(ms_maml + 1e-12), '-', 'LineWidth', 1.8, 'Color', [0.9290 0.6940 0.1250]);
plot(t, 10*log10(ms_aux_original + 1e-12), ':', 'LineWidth', 1.6, 'Color', [0.4940 0.1840 0.5560]);
plot(t, 10*log10(ms_aux_maml + 1e-12), ':', 'LineWidth', 1.6, 'Color', [0.20 0.35 0.85]);
xline(5, 'k:', 'Phase 1 (S_1, P_1)', ...
    'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
xline(t_switch1, 'k:', 'Phase 2 (switch to S_2, P_2)', ...
    'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
xline(t_switch2, 'k:', 'Phase 3 (switch to S_3, P_3)', ...
    'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
xlabel('Time (s)');
ylabel('Mean Squared Error (dB)');
title('Online modeling FxLMS with auxiliary-noise power (switches at 60 s & 120 s)');
legend({'ANC off', 'Error (Original)', 'Error (Proposed MAML)', ...
        'Aux noise (Original)', 'Aux noise (Proposed MAML)'}, ...
        'Location', 'best');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlim([0 sum(phase_sec)]);
ylim([-70 -10]);

save_figure(fig, fullfile(out_dir, 'Fig2_online_switch.png'));

function [e, st] = run_osmp_phase(ref_noisy, d_noisy, S_true, ...
                                  w0, s0, h0, ...
                                  mu_w, mu_s, mu_h, Tw, Ts, Th, ...
                                  alpha_power, c_aux, delay)
    Lw = numel(w0);
    Ls = numel(s0);
    Lh = numel(h0);
    n_samp = numel(ref_noisy);

    e = zeros(n_samp, 1);
    aux_mic = zeros(n_samp, 1);
    xw = zeros(Lw, 1);
    xh = zeros(Lh, 1);
    xs = zeros(Ls + Lw, 1);
    u_hist = zeros(Ls, 1);
    v_hist = zeros(Ls, 1);

    w = w0(:);
    s_hat = s0(:);
    h = h0(:);
    Px = 1;
    Pe1_power = 0;

    for n = 1:n_samp
        xw = [ref_noisy(n); xw(1:end-1)];
        if n < delay + 1
            xh = [0; xh(1:end-1)];
        else
            xh = [ref_noisy(n-delay); xh(1:end-1)];
        end
        xs = [xs(2:end); ref_noisy(n)];

        u_hist = [w.' * xw; u_hist(1:end-1)];

        probe = randn();
        if Pe1_power < Px
            v = c_aux * probe * sqrt(Pe1_power + 1e-12);
        else
            v = c_aux * probe * sqrt(Px + 1e-12);
        end
        v_hist = [v; v_hist(1:end-1)];

        y_control = u_hist.' * S_true;
        y_aux = v_hist.' * S_true;
        e(n) = d_noisy(n) - y_control + y_aux;
        aux_mic(n) = y_aux;

        v_hat = v_hist.' * s_hat;
        e_clean = e(n) - v_hat;
        z_hat = xh.' * h;
        e_h = e_clean - z_hat;
        e_s = e(n) - z_hat - v_hat;

        s_new = s_hat + mu_s * (v_hist * e_s);
        if norm(s_new) <= Ts
            s_hat = s_new;
        else
            s_hat = Ts * s_new / (norm(s_new) + 1e-12);
        end

        fx_buf = filter(s_hat, 1, xs);
        x_filtered = fx_buf(end:-1:end-Lw+1);
        w_new = w + mu_w * x_filtered * e_clean;
        if norm(w_new) <= Tw
            w = w_new;
        else
            w = Tw * w_new / (norm(w_new) + 1e-12);
        end

        h_new = h + mu_h * xh * e_h;
        if norm(h_new) <= Th
            h = h_new;
        else
            h = Th * h_new / (norm(h_new) + 1e-12);
        end

        Px = alpha_power * Px + (1 - alpha_power) * ref_noisy(n)^2;
        Pe1_power = alpha_power * Pe1_power + (1 - alpha_power) * e_clean^2;
    end

    st.w = w;
    st.s = s_hat;
    st.h = h;
    st.aux = aux_mic;
end

function [P, S] = pick_pair(meta, Hpri, Hsec, row, fs, Fs0)
    P = Hpri(:, meta.prim_col(row));
    S = Hsec(:, meta.sec_col(row));
    if fs ~= Fs0
        P = resample(P, fs, Fs0);
        S = resample(S, fs, Fs0);
    end
end

function xL = pad_or_trim(x, L)
    x = x(:);
    if numel(x) >= L
        xL = x(1:L);
    else
        xL = [x; zeros(L - numel(x), 1)];
    end
end

function y = add_awgn_measured(x, snr_db)
    p_signal = mean(x(:).^2);
    p_noise = p_signal / (10^(snr_db/10));
    y = x + sqrt(p_noise) * randn(size(x));
end

function save_figure(fig, filename)
    try
        exportgraphics(fig, filename, 'Resolution', 300);
    catch
        saveas(fig, filename);
    end
    fprintf('Saved %s\n', filename);
end
