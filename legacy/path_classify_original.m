%% === RWTH（out_pairs）次级路径汇总 + LSD分组（ALL ears, no delay alignment） ===
clear; clc;

%=== 可调参数 ===%
PAIRS_FILE = 'ANC_PriSec_pairs.mat';   % 数据文件
Nfft       = 4096;                     % 频点数
F_BAND     = [100 2000];               % LSD计算频带
AMP_MODE   = 'global';                 % 'absolute' | 'global' | 'medianrel'
DELAY_ALIGN_BY_PEAK = false;           % ★ 不对齐时延（按你要求）
GRAY = [0.7 0.7 0.7];                  % 灰线颜色
M_LSD_PTS  = 64;                       % LSD对数频轴采样点数
eps0 = 1e-12;

%=== 读取数据 ===%
S     = load(PAIRS_FILE, 'out_pairs');
Hsec  = S.out_pairs.H_secondary;      % [Nsamp x Ncols]
meta  = S.out_pairs.meta_pairs;       % table：含 sec_col/side/fit 等
Fs0   = S.out_pairs.Fs;

%=== 取全部有效次级路径索引（不过滤左右耳/佩戴） ===%
sec_col_all = meta.sec_col;
valid = isfinite(sec_col_all) & sec_col_all>=1 & sec_col_all<=size(Hsec,2);
sec_col = unique(sec_col_all(valid));                 % 去重
assert(~isempty(sec_col), '没有发现有效的次级路径索引 sec_col。');

%=== 频响计算（不做时延对齐） ===%
nSel  = numel(sec_col);
magAll = []; phAll = []; f_ref = [];
for k = 1:nSel
    c = sec_col(k);
    h = double(Hsec(:, c));                          % [Nsamp x 1]
    [H, f]  = freqz(h, 1, Nfft, Fs0);               % 原始H(f)（不补偿时延）
    magAll(:,k) = abs(H);                            %#ok<AGROW>
    phAll(:,k)  = unwrap(angle(H))*180/pi;           %#ok<AGROW>
    if isempty(f_ref), f_ref = f; end
end
maskBand = (f_ref>=F_BAND(1) & f_ref<=F_BAND(2));

%=== 中位线 & 统一参考（用于可视化） ===%
magMed   = median(magAll, 2, 'omitnan');
phMed    = median(phAll,   2, 'omitnan');
globalRef = max(magAll(maskBand,:), [], 'all');

%% === 计算 LSD 距离矩阵（不对齐时延） ===
f_log = logspace(log10(F_BAND(1)), log10(F_BAND(2)), M_LSD_PTS).';
XdB = zeros(nSel, numel(f_log));
for i = 1:nSel
    Mag_dB = 20*log10(magAll(:,i)+eps0);
    XdB(i,:) = interp1(f_ref, Mag_dB, f_log, 'linear', 'extrap');
end
% LSD = dB向量RMSE
D = zeros(nSel);
for i = 1:nSel
    diff = XdB(i,:) - XdB;            % [nSel, M_LSD_PTS]
    D(i,:) = sqrt(mean(diff.^2, 2));
end

%% === 选三条：A组（差异最大, maximin）& B组（差异最小, min-diameter） ===
bestA_score = -inf;  A_idx = [1 2 3]; A_dists = [NaN NaN NaN];
bestB_score =  inf;  B_idx = [1 2 3]; B_dists = [NaN NaN NaN];

for i = 1:nSel
    for j = i+1:nSel
        for k = j+1:nSel
            d1 = D(i,j); d2 = D(i,k); d3 = D(j,k);
            % A组指标：min pairwise LSD（越大越好）
            sA = min([d1 d2 d3]);
            if sA > bestA_score
                bestA_score = sA; A_idx = [i j k]; A_dists = [d1 d2 d3];
            end
            % B组指标：max pairwise LSD（越小越好）= 集合直径
            sB = max([d1 d2 d3]);
            if sB < bestB_score
                bestB_score = sB; B_idx = [i j k]; B_dists = [d1 d2 d3];
            end
        end
    end
end

% 打印结果（评价指标）
fprintf('\n=== LSD-based triad selection (no delay alignment) ===\n');
fprintf('A组(差异最大, maximin): idx=%s (sec_col=%s)\n', mat2str(A_idx), mat2str(sec_col(A_idx)));
fprintf('  pairwise LSDs (dB) = [%.3f %.3f %.3f], 指标(min)=%.3f dB\n', A_dists(1),A_dists(2),A_dists(3), bestA_score);
fprintf('B组(差异最小, min-diameter): idx=%s (sec_col=%s)\n', mat2str(B_idx), mat2str(sec_col(B_idx)));
fprintf('  pairwise LSDs (dB) = [%.3f %.3f %.3f], 指标(max)=%.3f dB\n', B_dists(1),B_dists(2),B_dists(3), bestB_score);

%% === 画幅度/相位 + 高亮两组 ===
figure('Name','RWTH - Secondary (InnerDriver -> Eardrum): ALL ears');
tiledlayout(2,1,'Padding','compact','TileSpacing','tight');
ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on');
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');

switch lower(AMP_MODE)
    case 'absolute'
        for k = 1:nSel
            plot(ax1, f_ref, 20*log10(magAll(:,k)+eps0), 'Color',GRAY, 'LineWidth',0.8);
        end
        medLine = 20*log10(magMed+eps0);
        plot(ax1, f_ref, medLine, 'k', 'LineWidth', 2.0);
        title(ax1, 'Magnitude (dB, absolute) — ALL ears');
    case 'global'
        ref = globalRef + eps0;
        for k = 1:nSel
            plot(ax1, f_ref, 20*log10(magAll(:,k)./ref + eps0), 'Color',GRAY, 'LineWidth',0.8);
        end
        medLine = 20*log10(magMed./ref + eps0);
        plot(ax1, f_ref, medLine, 'k', 'LineWidth', 2.0);
        title(ax1, sprintf('Magnitude (dB re global max in [%d,%d] Hz) — ALL ears', F_BAND(1), F_BAND(2)));
    case 'medianrel'
        med_dB = 20*log10(magMed+eps0);
        for k = 1:nSel
            plot(ax1, f_ref, 20*log10(magAll(:,k)+eps0)-med_dB, 'Color',GRAY, 'LineWidth',0.8);
        end
        medLine = zeros(size(f_ref));
        plot(ax1, f_ref, medLine, 'k', 'LineWidth', 2.0);
        title(ax1, 'ΔMagnitude vs median (dB) — ALL ears');
    otherwise
        error('AMP_MODE 必须是 ''absolute'' | ''global'' | ''medianrel''。');
end
ylabel(ax1,'Mag (dB)');

% 相位（灰线+中位）
for k = 1:nSel
    plot(ax2, f_ref, phAll(:,k), 'Color',GRAY, 'LineWidth',0.8);
end
plot(ax2, f_ref, phMed, 'k', 'LineWidth', 2.0);
title(ax2, 'Phase (deg, no delay alignment) — ALL ears');
xlabel(ax2,'Frequency (Hz)'); ylabel(ax2,'Phase (deg)');
xlim(ax1, [0, Fs0/2]); xlim(ax2, [0, Fs0/2]);

% 在幅度图上高亮 A/B 两组
colA = [0.85 0.20 0.20];  % 红
colB = [0.20 0.45 0.85];  % 蓝
for kk = 1:numel(A_idx)
    i = A_idx(kk);
    switch lower(AMP_MODE)
        case 'absolute'
            y = 20*log10(magAll(:,i)+eps0);
        case 'global'
            y = 20*log10(magAll(:,i)./(globalRef+eps0) + eps0);
        case 'medianrel'
            y = 20*log10(magAll(:,i)+eps0) - 20*log10(magMed+eps0);
    end
    plot(ax1, f_ref, y, 'Color',colA, 'LineWidth',2.2);
    plot(ax2, f_ref, phAll(:,i), 'Color',colA, 'LineWidth',2.2);
end
for kk = 1:numel(B_idx)
    i = B_idx(kk);
    switch lower(AMP_MODE)
        case 'absolute'
            y = 20*log10(magAll(:,i)+eps0);
        case 'global'
            y = 20*log10(magAll(:,i)./(globalRef+eps0) + eps0);
        case 'medianrel'
            y = 20*log10(magAll(:,i)+eps0) - 20*log10(magMed+eps0);
    end
    plot(ax1, f_ref, y, 'Color',colB, 'LineWidth',2.2, 'LineStyle','--');
    plot(ax2, f_ref, phAll(:,i), 'Color',colB, 'LineWidth',2.2, 'LineStyle','--');
end

% 目标频带阴影
yl = ylim(ax1);
patch(ax1, [F_BAND(1) F_BAND(2) F_BAND(2) F_BAND(1)], [yl(1) yl(1) yl(2) yl(2)], ...
      [0.9 0.95 1.0], 'EdgeColor','none', 'FaceAlpha', 0.15, 'HitTest','off');
uistack(ax1.Children(1),'bottom');

% 简洁图例
legend(ax1, {sprintf('median (N=%d)', nSel), 'A组(差异最大)', 'B组(差异最小)'}, 'Location','southwest');


%% === RWTH（out_pairs）初级路径：LSD分组（ALL ears, no delay alignment） ===

% 读取 primary
Hpri = S.out_pairs.H_primary;

% 取全部有效初级路径索引（不过滤左右耳/佩戴）
assert(ismember('prim_col', meta.Properties.VariableNames), ...
       'meta_pairs 中未找到字段 pri_col。');
pri_col_all = meta.prim_col;
valid_p = isfinite(pri_col_all) & pri_col_all>=1 & pri_col_all<=size(Hpri,2);
pri_col = unique(pri_col_all(valid_p));
assert(numel(pri_col) >= 3, '有效的初级路径不足 3 条。');

% 频响计算（不做时延对齐）
nP = numel(pri_col);
magAll_p = []; phAll_p = []; f_ref_p = [];
for k = 1:nP
    c = pri_col(k);
    h = double(Hpri(:, c));
    [H, f] = freqz(h, 1, Nfft, Fs0);      % 原始H(f)（不补偿时延）
    magAll_p(:,k) = abs(H);               %#ok<AGROW>
    phAll_p(:,k)  = unwrap(angle(H))*180/pi; %#ok<AGROW>
    if isempty(f_ref_p), f_ref_p = f; end
end
maskBand_p = (f_ref_p>=F_BAND(1) & f_ref_p<=F_BAND(2));

% 中位线 & 统一参考（用于可视化）
magMed_p   = median(magAll_p, 2, 'omitnan');
phMed_p    = median(phAll_p,   2, 'omitnan');
globalRef_p = max(magAll_p(maskBand_p,:), [], 'all');

% 计算 LSD 距离矩阵（不对齐时延）
f_log = logspace(log10(F_BAND(1)), log10(F_BAND(2)), M_LSD_PTS).';
XdB_p = zeros(nP, numel(f_log));
for i = 1:nP
    Mag_dB = 20*log10(magAll_p(:,i)+eps0);
    XdB_p(i,:) = interp1(f_ref_p, Mag_dB, f_log, 'linear', 'extrap');
end
D_p = zeros(nP);
for i = 1:nP
    diff = XdB_p(i,:) - XdB_p;           % [nP, M_LSD_PTS]
    D_p(i,:) = sqrt(mean(diff.^2, 2));
end

% 选三条：A组（差异最大, maximin）& B组（差异最小, min-diameter）
bestA_p = -inf;  A_idx_p = [1 2 3]; A_dists_p = [NaN NaN NaN];
bestB_p =  inf;  B_idx_p = [1 2 3]; B_dists_p = [NaN NaN NaN];
for i = 1:nP
    for j = i+1:nP
        for k = j+1:nP
            d1 = D_p(i,j); d2 = D_p(i,k); d3 = D_p(j,k);
            sA = min([d1 d2 d3]);
            if sA > bestA_p
                bestA_p = sA; A_idx_p = [i j k]; A_dists_p = [d1 d2 d3];
            end
            sB = max([d1 d2 d3]);
            if sB < bestB_p
                bestB_p = sB; B_idx_p = [i j k]; B_dists_p = [d1 d2 d3];
            end
        end
    end
end

% 打印结果（评价指标）
fprintf('\n=== PRIMARY: LSD-based triad (no delay alignment) ===\n');
fprintf('A组(差异最大, maximin): idx=%s (pri_col=%s)\n', ...
        mat2str(A_idx_p), mat2str(pri_col(A_idx_p)));
fprintf('  pairwise LSDs (dB) = [%.3f %.3f %.3f], 指标(min)=%.3f dB\n', ...
        A_dists_p(1),A_dists_p(2),A_dists_p(3), bestA_p);
fprintf('B组(差异最小, min-diameter): idx=%s (pri_col=%s)\n', ...
        mat2str(B_idx_p), mat2str(pri_col(B_idx_p)));
fprintf('  pairwise LSDs (dB) = [%.3f %.3f %.3f], 指标(max)=%.3f dB\n', ...
        B_dists_p(1),B_dists_p(2),B_dists_p(3), bestB_p);

% 画幅度/相位 + 高亮两组
figure('Name','RWTH - Primary (noise -> error): ALL ears');
tiledlayout(2,1,'Padding','compact','TileSpacing','tight');
ax1p = nexttile; hold(ax1p,'on'); grid(ax1p,'on');
ax2p = nexttile; hold(ax2p,'on'); grid(ax2p,'on');

switch lower(AMP_MODE)
    case 'absolute'
        for k = 1:nP
            plot(ax1p, f_ref_p, 20*log10(magAll_p(:,k)+eps0), 'Color',GRAY, 'LineWidth',0.8);
        end
        medLine_p = 20*log10(magMed_p+eps0);
        plot(ax1p, f_ref_p, medLine_p, 'k', 'LineWidth', 2.0);
        title(ax1p, 'Magnitude (dB, absolute) — ALL ears');
    case 'global'
        refp = globalRef_p + eps0;
        for k = 1:nP
            plot(ax1p, f_ref_p, 20*log10(magAll_p(:,k)./refp + eps0), 'Color',GRAY, 'LineWidth',0.8);
        end
        medLine_p = 20*log10(magMed_p./refp + eps0);
        plot(ax1p, f_ref_p, medLine_p, 'k', 'LineWidth', 2.0);
        title(ax1p, sprintf('Magnitude (dB re global max in [%d,%d] Hz) — ALL ears', F_BAND(1), F_BAND(2)));
    case 'medianrel'
        med_dB_p = 20*log10(magMed_p+eps0);
        for k = 1:nP
            plot(ax1p, f_ref_p, 20*log10(magAll_p(:,k)+eps0)-med_dB_p, 'Color',GRAY, 'LineWidth',0.8);
        end
        plot(ax1p, f_ref_p, zeros(size(f_ref_p)), 'k', 'LineWidth', 2.0);
        title(ax1p, 'ΔMagnitude vs median (dB) — ALL ears');
end
ylabel(ax1p,'Mag (dB)');

for k = 1:nP
    plot(ax2p, f_ref_p, phAll_p(:,k), 'Color',GRAY, 'LineWidth',0.8);
end
plot(ax2p, f_ref_p, phMed_p, 'k', 'LineWidth', 2.0);
title(ax2p, 'Phase (deg, no delay alignment) — ALL ears');
xlabel(ax2p,'Frequency (Hz)'); ylabel(ax2p,'Phase (deg)');
xlim(ax1p, [0, Fs0/2]); xlim(ax2p, [0, Fs0/2]);

colA = [0.85 0.20 0.20];  % 红
colB = [0.20 0.45 0.85];  % 蓝
for kk = 1:numel(A_idx_p)
    i = A_idx_p(kk);
    switch lower(AMP_MODE)
        case 'absolute'
            y = 20*log10(magAll_p(:,i)+eps0);
        case 'global'
            y = 20*log10(magAll_p(:,i)./(globalRef_p+eps0) + eps0);
        case 'medianrel'
            y = 20*log10(magAll_p(:,i)+eps0) - 20*log10(magMed_p+eps0);
    end
    plot(ax1p, f_ref_p, y, 'Color',colA, 'LineWidth',2.2);
    plot(ax2p, f_ref_p, phAll_p(:,i), 'Color',colA, 'LineWidth',2.2);
end
for kk = 1:numel(B_idx_p)
    i = B_idx_p(kk);
    switch lower(AMP_MODE)
        case 'absolute'
            y = 20*log10(magAll_p(:,i)+eps0);
        case 'global'
            y = 20*log10(magAll_p(:,i)./(globalRef_p+eps0) + eps0);
        case 'medianrel'
            y = 20*log10(magAll_p(:,i)+eps0) - 20*log10(magMed_p+eps0);
    end
    plot(ax1p, f_ref_p, y, 'Color',colB, 'LineWidth',2.2, 'LineStyle','--');
    plot(ax2p, f_ref_p, phAll_p(:,i), 'Color',colB, 'LineWidth',2.2, 'LineStyle','--');
end

yl = ylim(ax1p);
patch(ax1p, [F_BAND(1) F_BAND(2) F_BAND(2) F_BAND(1)], [yl(1) yl(1) yl(2) yl(2)], ...
      [0.9 0.95 1.0], 'EdgeColor','none', 'FaceAlpha', 0.15, 'HitTest','off');
uistack(ax1p.Children(1),'bottom');
legend(ax1p, {sprintf('median (N=%d)', nP), 'A组(差异最大)', 'B组(差异最小)'}, 'Location','southwest');


%% === Row-aligned triads: select on Secondary ↔ report Primary, and vice versa ===

Hpri = S.out_pairs.H_primary;
Hsec = S.out_pairs.H_secondary;

% —— 以 meta 的“行”为单位，确保一一对应（每行有 pri_col & sec_col）——
assert(all(ismember({'prim_col','sec_col'}, meta.Properties.VariableNames)), ...
       'meta_pairs 需要包含 pri_col 与 sec_col 字段。');

rows_valid = isfinite(meta.prim_col) & meta.prim_col>=1 & meta.prim_col<=size(Hpri,2) & ...
             isfinite(meta.sec_col) & meta.sec_col>=1 & meta.sec_col<=size(Hsec,2);
rows = find(rows_valid);
nr   = numel(rows);
assert(nr>=3, '有效行不足 3 条。');

% —— 计算每一行的 Primary / Secondary 幅度谱（不做时延对齐）——
magP = zeros(Nfft, nr);
magS = zeros(Nfft, nr);
for t = 1:nr
    r = rows(t);
    [HP, f_ref2] = freqz(double(Hpri(:, meta.prim_col(r))), 1, Nfft, Fs0);
    [HS, ~     ] = freqz(double(Hsec(:, meta.sec_col(r))), 1, Nfft, Fs0);
    magP(:,t) = abs(HP);
    magS(:,t) = abs(HS);
end

% —— 两侧各自的 LSD 距离矩阵（同一组行索引上可互查）——
D_pri = lsd_dist_from_mag(magP, f_ref2, F_BAND, M_LSD_PTS, eps0);
D_sec = lsd_dist_from_mag(magS, f_ref2, F_BAND, M_LSD_PTS, eps0);

% ================== ① 按“次级路径”选，然后报告初级路径的 LSD ==================
[A_idx_s, A_score_s, A_trip_s] = triad_maximin(D_sec);
[B_idx_s, B_score_s, B_trip_s] = triad_min_diameter(D_sec);

print_triad('按次级路径选择  A组(差异最大, maximin)', ...
    rows, meta, D_sec, D_pri, A_idx_s, A_score_s, A_trip_s, true);
print_triad('按次级路径选择  B组(差异最小, min-diameter)', ...
    rows, meta, D_sec, D_pri, B_idx_s, B_score_s, B_trip_s, false);

% ================== ② 按“初级路径”选，然后报告次级路径的 LSD ==================
[A_idx_p, A_score_p, A_trip_p] = triad_maximin(D_pri);
[B_idx_p, B_score_p, B_trip_p] = triad_min_diameter(D_pri);

print_triad('按初级路径选择  A组(差异最大, maximin)', ...
    rows, meta, D_pri, D_sec, A_idx_p, A_score_p, A_trip_p, true);
print_triad('按初级路径选择  B组(差异最小, min-diameter)', ...
    rows, meta, D_pri, D_sec, B_idx_p, B_score_p, B_trip_p, false);

%% ======== 本段用到的极简函数 ========
function D = lsd_dist_from_mag(magAll, f_ref, F_BAND, M, eps0)
    % 将每条曲线映射到对数频轴上的 dB 向量，再做 RMSE
    f_log = logspace(log10(F_BAND(1)), log10(F_BAND(2)), M).';
    N = size(magAll,2);
    XdB = zeros(N, numel(f_log));
    for i = 1:N
        XdB(i,:) = interp1(f_ref, 20*log10(magAll(:,i)+eps0), f_log, 'linear','extrap');
    end
    D = zeros(N);
    for i = 1:N
        diff = XdB(i,:) - XdB;                  % [N, M]
        D(i,:) = sqrt(mean(diff.^2, 2));        % LSD(dB)
    end
end

function [idx, score, dists] = triad_maximin(D)
    % 选三条，使 min(pairwise D) 最大
    N = size(D,1); score = -inf; idx=[1 2 3]; dists=[NaN NaN NaN];
    for i=1:N, for j=i+1:N, for k=j+1:N
        d = [D(i,j) D(i,k) D(j,k)]; s = min(d);
        if s > score, score = s; idx = [i j k]; dists = d; end
    end, end, end
end

function [idx, score, dists] = triad_min_diameter(D)
    % 选三条，使 max(pairwise D) 最小（直径最小）
    N = size(D,1); score = inf; idx=[1 2 3]; dists=[NaN NaN NaN];
    for i=1:N, for j=i+1:N, for k=j+1:N
        d = [D(i,j) D(i,k) D(j,k)]; s = max(d);
        if s < score, score = s; idx = [i j k]; dists = d; end
    end, end, end
end

function print_triad(title_str, rows, meta, D_main, D_cross, idx, score, dists, isA)
    % D_main: 用来“选”的那一侧；D_cross: 另一侧（用于同时报告）
    i=idx(1); j=idx(2); k=idx(3);
    rsel = rows(idx);
    % 对应列号
    pri_cols = meta.prim_col(rsel);
    sec_cols = meta.sec_col(rsel);

    % 另一侧的三对 LSD
    dc = [D_cross(i,j) D_cross(i,k) D_cross(j,k)];
    mean_main  = mean(dists); min_main = min(dists); max_main = max(dists);
    mean_cross = mean(dc);    min_cross = min(dc);   max_cross = max(dc);

    fprintf('\n=== %s ===\n', title_str);
    fprintf('rows = %s | pri_col = %s | sec_col = %s\n', ...
            mat2str(rsel), mat2str(pri_cols.'), mat2str(sec_cols.'));
    if isA
        fprintf('选用侧 LSD(dB): [% .3f % .3f % .3f] | min=%.3f ←主指标 | mean=%.3f | max=%.3f\n', ...
                dists(1), dists(2), dists(3), min_main, mean_main, max_main);
    else
        fprintf('选用侧 LSD(dB): [% .3f % .3f % .3f] | max=%.3f ←主指标 | mean=%.3f | min=%.3f\n', ...
                dists(1), dists(2), dists(3), max_main, mean_main, min_main);
    end
    fprintf('另一侧 LSD(dB):   [% .3f % .3f % .3f] | min=%.3f | mean=%.3f | max=%.3f\n', ...
            dc(1), dc(2), dc(3), min_cross, mean_cross, max_cross);
end



%% === ③ 以 “mean(LSD_primary)+mean(LSD_secondary)” 为准则的三元组排序（升序） ===
TOPK = 8000;                                   % 想看前多少个
nr = size(D_pri,1);
ncomb = nchoosek(nr,3);

score  = zeros(ncomb,1);
meanP  = zeros(ncomb,1);
meanS  = zeros(ncomb,1);
trip   = zeros(ncomb,3);

t = 0;
for i = 1:nr-2
    for j = i+1:nr-1
        for k = j+1:nr
            t = t + 1;
            dp = [D_pri(i,j), D_pri(i,k), D_pri(j,k)];
            ds = [D_sec(i,j), D_sec(i,k), D_sec(j,k)];
            meanP(t) = mean(dp);
            meanS(t) = mean(ds);
            score(t) = meanP(t) + meanS(t);
            trip(t,:) = [i j k];
        end
    end
end

[score_sorted, ord] = sort(score, 'ascend');
K = min(TOPK, numel(ord));

fprintf('\n=== 按  mean(LSD_pri)+mean(LSD_sec)  升序的前 %d 个三元组 ===\n', K);
fprintf('   rank | rows           | pri_col        | sec_col        | meanP  | meanS  | sum\n');
fprintf('--------+----------------+----------------+----------------+--------+--------+--------\n');
for r = 1:K
    idx = trip(ord(r),:);
    rr  = rows(idx).';                               % 3×1 -> 1×3
    pc  = meta.prim_col(rr).';
    sc  = meta.sec_col(rr).';
    fprintf(' %6d | [%2d %2d %2d] | [%2d %2d %2d] | [%2d %2d %2d] | %6.3f | %6.3f | %6.3f\n', ...
        r, rr(1), rr(2), rr(3), pc(1), pc(2), pc(3), sc(1), sc(2), sc(3), ...
        meanP(ord(r)), meanS(ord(r)), score_sorted(r));
end

%% === 补充：总共有多少三元组 & 打印分数最大的 TOPK 组（降序） ===
fprintf('\n总组合数（从 %d 条有效行中取 3 条） = %d 组\n', nr, ncomb);

TOPK_HI = TOPK;                        % 想看多少个最大的；也可单独设
ord_desc = flipud(ord);                % 把升序索引翻转为降序
Khi = min(TOPK_HI, numel(ord_desc));

fprintf('\n=== 按  mean(LSD_pri)+mean(LSD_sec)  降序的前 %d 个三元组（“最大”）===\n', Khi);
fprintf('   rank | rows           | pri_col        | sec_col        | meanP  | meanS  | sum\n');
fprintf('--------+----------------+----------------+----------------+--------+--------+--------\n');
for rnk = 1:Khi
    idx = trip(ord_desc(rnk),:);
    rr  = rows(idx).';
    pc  = meta.prim_col(rr).';
    sc  = meta.sec_col(rr).';
    fprintf(' %6d | [%2d %2d %2d] | [%2d %2d %2d] | [%2d %2d %2d] | %6.3f | %6.3f | %6.3f\n', ...
        rnk, rr(1), rr(2), rr(3), pc(1), pc(2), pc(3), sc(1), sc(2), sc(3), ...
        meanP(ord_desc(rnk)), meanS(ord_desc(rnk)), score_sorted(end-rnk+1));
end

%% === 根据三条"行号"画图：先灰色全部，再叠加A/B两组（不画中位线） ===
% 把你选出来的行号填这里；B组不需要可设 []
ROWS_A = [32 37 43];   % A组三条 meta 行号
ROWS_B = [20 26 31];   % B组三条 meta 行号

% —— 调用新的绘图函数，一次性画Primary和Secondary ——
plot_primary_secondary_responses(S.out_pairs.H_primary, S.out_pairs.H_secondary, ...
    'prim_col', 'sec_col', ROWS_A, ROWS_B, ...
    Fs0, Nfft, F_BAND, AMP_MODE, eps0, GRAY, meta);

% ====================== 主绘图函数 ======================
function plot_primary_secondary_responses(H_primary, H_secondary, prim_col_field, sec_col_field, rowsA, rowsB, Fs0, Nfft, F_BAND, AMP_MODE, eps0, GRAY, meta)
    
    % --- 计算Primary和Secondary的频响数据 ---
    [magAll_prim, phAll_prim, f_ref, globalRef_prim, med_dB_prim] = compute_frequency_responses(H_primary, prim_col_field, meta, Fs0, Nfft, F_BAND, eps0);
    [magAll_sec, phAll_sec, ~, globalRef_sec, med_dB_sec] = compute_frequency_responses(H_secondary, sec_col_field, meta, Fs0, Nfft, F_BAND, eps0);
    
    % --- 第一张图：Primary和Secondary幅度响应 ---
    figure('Name','RWTH - Primary & Secondary — Magnitude Responses');
    tiledlayout(2,1,'Padding','compact','TileSpacing','tight');
    
    % Primary幅度响应
    ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on');
    plot_magnitude_responses(ax1, f_ref, magAll_prim, AMP_MODE, globalRef_prim, med_dB_prim, eps0, F_BAND, GRAY);
    if ~isempty(rowsA), overlay_magnitude_group(ax1, H_primary, meta.(prim_col_field)(rowsA), ...
            Fs0, Nfft, f_ref, AMP_MODE, globalRef_prim, med_dB_prim, eps0, [0.85 0.20 0.20], '-'); end
    if ~isempty(rowsB), overlay_magnitude_group(ax1, H_primary, meta.(prim_col_field)(rowsB), ...
            Fs0, Nfft, f_ref, AMP_MODE, globalRef_prim, med_dB_prim, eps0, [0.20 0.45 0.85], '--'); end
    add_frequency_band_shading(ax1, F_BAND);
    add_legend(ax1, rowsA, rowsB);
    title(ax1, 'Primary paths');
    xlim(ax1, [0, 2000]);
    ylim(ax1, [-50, 20]);  % 设置Primary路径的y轴范围
    
    % Secondary幅度响应
    ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');
    plot_magnitude_responses(ax2, f_ref, magAll_sec, AMP_MODE, globalRef_sec, med_dB_sec, eps0, F_BAND, GRAY);
    if ~isempty(rowsA), overlay_magnitude_group(ax2, H_secondary, meta.(sec_col_field)(rowsA), ...
            Fs0, Nfft, f_ref, AMP_MODE, globalRef_sec, med_dB_sec, eps0, [0.85 0.20 0.20], '-'); end
    if ~isempty(rowsB), overlay_magnitude_group(ax2, H_secondary, meta.(sec_col_field)(rowsB), ...
            Fs0, Nfft, f_ref, AMP_MODE, globalRef_sec, med_dB_sec, eps0, [0.20 0.45 0.85], '--'); end
    add_frequency_band_shading(ax2, F_BAND);
    add_legend(ax2, rowsA, rowsB);
    title(ax2, 'Secondary paths');
    xlabel(ax2,'Frequency (Hz)');
    xlim(ax2, [0, 2000]);

    % --- 第二张图：Primary和Secondary相位响应 ---
    figure('Name','RWTH - Primary & Secondary — Phase Responses');
    tiledlayout(2,1,'Padding','compact','TileSpacing','tight');
    
    % Primary相位响应
    ax3 = nexttile; hold(ax3,'on'); grid(ax3,'on');
    plot_phase_responses(ax3, f_ref, phAll_prim, GRAY);
    if ~isempty(rowsA), overlay_phase_group(ax3, H_primary, meta.(prim_col_field)(rowsA), ...
            Fs0, Nfft, f_ref, [0.85 0.20 0.20], '-'); end
    if ~isempty(rowsB), overlay_phase_group(ax3, H_primary, meta.(prim_col_field)(rowsB), ...
            Fs0, Nfft, f_ref, [0.20 0.45 0.85], '--'); end
    add_legend(ax3, rowsA, rowsB);
    title(ax3, 'Primary path');
    ylabel(ax3,'Phase (deg)');
    xlim(ax3, [0, 2000]);
    
    % Secondary相位响应
    ax4 = nexttile; hold(ax4,'on'); grid(ax4,'on');
    plot_phase_responses(ax4, f_ref, phAll_sec, GRAY);
    if ~isempty(rowsA), overlay_phase_group(ax4, H_secondary, meta.(sec_col_field)(rowsA), ...
            Fs0, Nfft, f_ref, [0.85 0.20 0.20], '-'); end
    if ~isempty(rowsB), overlay_phase_group(ax4, H_secondary, meta.(sec_col_field)(rowsB), ...
            Fs0, Nfft, f_ref, [0.20 0.45 0.85], '--'); end
    add_legend(ax4, rowsA, rowsB);
    title(ax4, 'Secondary Paths - Phase');
    xlabel(ax4,'Frequency (Hz)'); ylabel(ax4,'Phase (deg)');
    xlim(ax4, [0, 2000]);
end

% ====================== 辅助函数 ======================

% --- 计算频响数据 ---
function [magAll, phAll, f_ref, globalRef, med_dB] = compute_frequency_responses(Hmat, col_field, meta, Fs0, Nfft, F_BAND, eps0)
    assert(ismember(col_field, meta.Properties.VariableNames), ...
        'meta_pairs 中未找到字段 %s', col_field);
    col_all = meta.(col_field);
    valid   = isfinite(col_all) & col_all>=1 & col_all<=size(Hmat,2);
    cols    = unique(col_all(valid));
    nAll    = numel(cols);
    assert(nAll>=1, '没有有效列可画。');

    magAll = zeros(Nfft, nAll);
    phAll  = zeros(Nfft, nAll);
    f_ref  = [];
    for i = 1:nAll
        h = double(Hmat(:, cols(i)));
        [H, f] = freqz(h, 1, Nfft, Fs0);
        magAll(:,i) = abs(H);
        phAll(:,i)  = unwrap(angle(H))*180/pi;
        if isempty(f_ref), f_ref = f; end
    end
    maskBand  = (f_ref>=F_BAND(1) & f_ref<=F_BAND(2));
    globalRef = max(magAll(maskBand,:), [], 'all');
    med_dB    = 20*log10(median(magAll,2,'omitnan')+eps0);
end

% --- 绘制幅度响应 ---
function plot_magnitude_responses(ax, f_ref, magAll, AMP_MODE, globalRef, med_dB, eps0, F_BAND, GRAY)
    nAll = size(magAll, 2);
    switch lower(AMP_MODE)
        case 'absolute'
            for i = 1:nAll
                plot(ax, f_ref, 20*log10(magAll(:,i)+eps0), 'Color',GRAY, 'LineWidth',0.6);
            end
            ylabel(ax,'Magnitude (dB)');
        case 'global'
            ref = globalRef + eps0;
            for i = 1:nAll
                plot(ax, f_ref, 20*log10(magAll(:,i)./ref + eps0), 'Color',GRAY, 'LineWidth',0.6);
            end
            ylabel(ax, sprintf('Magnitude (dB)', F_BAND(1), F_BAND(2)));
        case 'medianrel'
            for i = 1:nAll
                plot(ax, f_ref, 20*log10(magAll(:,i)+eps0) - med_dB, 'Color',GRAY, 'LineWidth',0.6);
            end
            ylabel(ax,'ΔMagnitude vs median (dB)');
        otherwise
            error('AMP_MODE 必须是 ''absolute'' | ''global'' | ''medianrel''。');
    end
end

% --- 绘制相位响应 ---
function plot_phase_responses(ax, f_ref, phAll, GRAY)
    nAll = size(phAll, 2);
    for i = 1:nAll
        plot(ax, f_ref, phAll(:,i), 'Color',GRAY, 'LineWidth',0.6);
    end
    ylabel(ax,'Phase (deg)');
end

% --- 叠加幅度组 ---
function overlay_magnitude_group(ax, Hmat, cols3, Fs0, Nfft, f_ref, AMP_MODE, globalRef, med_dB, eps0, colorRGB, style)
    cols3 = cols3(:).';
    for k = 1:numel(cols3)
        h = double(Hmat(:, cols3(k)));
        H = freqz(h, 1, Nfft, Fs0);
        mag = abs(H);
        switch lower(AMP_MODE)
            case 'absolute'
                y = 20*log10(mag+eps0);
            case 'global'
                y = 20*log10(mag./(globalRef+eps0) + eps0);
            case 'medianrel'
                y = 20*log10(mag+eps0) - med_dB;
        end
        plot(ax, f_ref, y, 'Color',colorRGB, 'LineWidth',1.8, 'LineStyle',style);
    end
end

% --- 叠加相位组 ---
function overlay_phase_group(ax, Hmat, cols3, Fs0, Nfft, f_ref, colorRGB, style)
    cols3 = cols3(:).';
    for k = 1:numel(cols3)
        h = double(Hmat(:, cols3(k)));
        H = freqz(h, 1, Nfft, Fs0);
        ph = unwrap(angle(H))*180/pi;
        plot(ax, f_ref, ph, 'Color',colorRGB, 'LineWidth',1.8, 'LineStyle',style);
    end
end

% --- 添加频带阴影 ---
function add_frequency_band_shading(ax, F_BAND)
    yl = ylim(ax);
    patch(ax, [F_BAND(1) F_BAND(2) F_BAND(2) F_BAND(1)], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.9 0.95 1.0], 'EdgeColor','none', 'FaceAlpha',0.15, 'HitTest','off');
    uistack(ax.Children(1),'bottom');
end

% --- 添加英文图例 ---
function add_legend(ax, rowsA, rowsB)
    % 创建图例条目
    lg_labels = {};
    lg_handles = [];
    
    % 灰色线（所有响应）
    lg_labels{end+1} = 'All responses';
    lg_handles(end+1) = plot(ax, NaN, NaN, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.6);
    
    % 红色线（Group A）
    if ~isempty(rowsA)
        lg_labels{end+1} = 'Group A';
        lg_handles(end+1) = plot(ax, NaN, NaN, 'Color', [0.85 0.20 0.20], 'LineWidth', 1.8, 'LineStyle', '-');
    end
    
    % 蓝色线（Group B）
    if ~isempty(rowsB)
        lg_labels{end+1} = 'Group B';
        lg_handles(end+1) = plot(ax, NaN, NaN, 'Color', [0.20 0.45 0.85], 'LineWidth', 1.8, 'LineStyle', '--');
    end
    
    % 创建图例
    legend(ax, lg_handles, lg_labels, 'Location', 'southwest');
end