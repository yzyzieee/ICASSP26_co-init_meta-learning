% close all;
clear; clc;

%% ============================ 配置参数 ============================
fs        = 16000;            % 目标采样率（保持与参考程序一致）
Len_N     = 512;              % 控制滤波器长度 = MAML Phi 长度
T_noise   = 3;                % 每次造噪声的时长（秒）
band_hz   = [100 1000];        % 与参考程序一致的带限
fir_ord   = 512;

Ls_S     = 256;     % 次级路径 FIR 长度
muS      = 0.05;     % 与W同风格可调
lamdaS   = 0.99;
epslonS  = 0.5;

% --- MAML 超参（保持参考程序的调用方式） ---
N_epcho   = 4096 * 20;       % 训练步数（按你原程序）
mu        = 0.1;
lamda     = 0.99;
epslon    = 0.5;

% ----------【选择数据】用 meta_pairs 的“行号”----------
PAIRS_FILE = 'ANC_PriSec_pairs.mat';   % 里面有 out_pairs 结构（你之前整理的）
% TRAIN_ROWS = [10 35 7];      % 训练任务池：随便放几个行号
% TEST_ROW   = 25;                       % 测试用哪一行
% BASE_ROW   = 10;                       % 单任务基线训练用哪一行
% TRAIN_ROWS = [32 37 43];      % 训练任务池：随便放几个行号
% TRAIN_ROWS = [20 26 31]; %min
% TRAIN_ROWS = [4 25 46]; %C
TRAIN_ROWS = [9 12 21]; %D
TEST_ROW   = 25;                       % 测试用哪一行
BASE_ROW   = 10;                       % 单任务基线训练用哪一行
% ---------------------------------------------------------------


% ====== 训练进度与可视化（NEW） ======
ema_rho    = 0.01;                          % EMA系数（仅用于画图）
S_pow_ema  = 0;  W_pow_ema = 0;             % 误差功率的EMA
S_pow_hist = zeros(N_epcho,1);
W_pow_hist = zeros(N_epcho,1);
print_every = max(1, floor(N_epcho/100));   % 每1%打印一次

use_waitbar = true; hwb = [];
if use_waitbar
    try hwb = waitbar(0,'Training Phase A/B ...'); catch, use_waitbar = false; end
end

rng(0);




%% ======================= 载入 RWTH 配对数据 =======================
S     = load(PAIRS_FILE,'out_pairs');
Fs0   = S.out_pairs.Fs;               % 48000
Hpri  = S.out_pairs.H_primary;        % [N x 46]
Hsec  = S.out_pairs.H_secondary;      % [N x 46]
meta  = S.out_pairs.meta_pairs;       % 表: person_id/side/fit/sec_col/prim_col

% 带限噪声生成器
bp = fir1(fir_ord, band_hz/(fs/2));
mk_noise = @(T) filter(bp,1, randn(round(T*fs),1));

% 小工具：按“行号”取成对(P,S)并重采样到 fs
get_pair = @(row) local_pick_pair(Hpri,Hsec,meta,row,fs,Fs0);

%% ====================== 准备训练样本（RWTH） ======================
% 关键替换：把参考程序里
%   x_ref = filter(P_ref,1,white); xprime = filter(S,1,x_ref); d = filter(P_err,1,white);
% 换成 RWTH 的
%   x  = 带限噪声;  d = filter(P,1,x);  xprime = filter(S,1,x)
% （没有 ref-mic 的情况下，用 x 自身做参考，符合你之前约定）

Fx_data = zeros(Len_N, N_epcho);
Di_data = zeros(Len_N, N_epcho);

train_rows = TRAIN_ROWS(:).';
assert(all(train_rows>=1 & train_rows<=height(meta)), '训练行号越界');

for jj = 1:N_epcho
    % 随机挑一个训练行 -> 成对(P,S)
    r   = train_rows(randi(numel(train_rows)));
    [P,S_pair,inf_r] = get_pair(r); %#ok<ASGLU>

    % 造带限噪声并通过路径
    x   = mk_noise(T_noise);
    d   = filter(P,1,x);          % 初级：噪声->误差
    xpf = filter(S_pair,1,x);     % 次级：参考->误差（用 x 作为参考）

    % 随机裁一段 Len_N（注意：外部不翻转，你的 MAML 类内部会 flipud）
    idx_cut = randi([Len_N, length(d)]);
    Di_data(:,jj) = d  (idx_cut-Len_N+1:idx_cut);
    Fx_data(:,jj) = xpf(idx_cut-Len_N+1:idx_cut);
end

%% ====================== 训练：Phase A (S) -> Phase B (W) ======================
aW = MAML_Nstep_forget(Len_N);       % 控制滤波器的元初始化 (已有类)
aS = MAML_Nstep_forget_S2(Ls_S);      % 次级路径的元初始化 (新增类)

ErS_train = zeros(N_epcho,1);        % 路径辨识误差（记录jj==1处的瞬时e，风格与W保持一致）
ErW_train = zeros(N_epcho,1);        % Error history

for jj = 1:N_epcho
    % -------- 随机抽任务行号（同你原来做法） --------
    r = train_rows(randi(numel(train_rows)));
    [P_true, S_true, ~] = get_pair(r);

    % ================= Phase A：次级路径辨识（固定3s段，2段音频：U和Y） =================
    % 生成辨识激励 U_id，并通过真实次级路径得到输出 Y_id
    U_id  = mk_noise(T_noise);            % excite
    Y_id  = filter(S_true,1,U_id);        % measured response

    % 裁取长度 Ls_S 的一段（翻转/对齐在类里做）
    cutS   = randi([Ls_S, length(Y_id)]);
    U_cut  = U_id(cutS-Ls_S+1 : cutS);
    Y_cut  = Y_id(cutS-Ls_S+1 : cutS);

    % 一步内环 + 外环元更新；得到该任务的 S_task
    [aS, S_task, ErS_train(jj)] = aS.MAML_initial_S(U_cut, Y_cut, muS, lamdaS, epslonS);

    % ================= Phase B：控制滤波器 (使用收敛后的 S_task) =================
    x   = mk_noise(T_noise);              % ANC 训练参考
    d   = filter(P_true,1,x);             % disturbance
    xp  = filter(S_task,1,x);             % 用 S_task 生成的 filtered-x

    % 裁取长度 Len_N 的一段
    cutW   = randi([Len_N, length(d)]);
    Di_cut = d (cutW-Len_N+1 : cutW);
    Fx_cut = xp(cutW-Len_N+1 : cutW);

    % 一步内环 + 外环元更新；W 的误差（记录风格与原来一致）
    [aW, ErW_train(jj)] = aW.MAML_initial(Fx_cut, Di_cut, mu, lamda, epslon);

        % ------ 记录/显示训练误差（NEW） ------
    % 用瞬时 e（你类里返回的 ErS_train(jj), ErW_train(jj)）做功率EMA，仅用于可视化
    S_pow_ema = (1-ema_rho)*S_pow_ema + ema_rho*(ErS_train(jj)^2);
    W_pow_ema = (1-ema_rho)*W_pow_ema + ema_rho*(ErW_train(jj)^2);
    S_pow_hist(jj) = S_pow_ema;
    W_pow_hist(jj) = W_pow_ema;
    
    if mod(jj, print_every)==0 || jj==N_epcho
        fprintf('[%5.1f%%%%]  S_EMA=%.2f dB | W_EMA=%.2f dB\n', ...
            100*jj/N_epcho, 10*log10(S_pow_ema+1e-12), 10*log10(W_pow_ema+1e-12));
        if use_waitbar && ~isempty(hwb)
            try waitbar(jj/N_epcho, hwb, sprintf('Epoch %d/%d', jj, N_epcho)); end
        end
    end
end

if use_waitbar && ~isempty(hwb), try close(hwb); end, end

figure;
plot(ErW_train); grid on;
xlabel('Epoch'); ylabel('Residual error');

figure;
plot(ErS_train); grid on;
xlabel('Epoch'); ylabel('Residual error');

Wc = aW.Phi;       % 控制滤波器的元初始化
Sc = aS.Psi;       % 次级路径的元初始化
save('Wc_Sc_MAML_RWTH_D_3.mat','Wc','Sc');

%% ========================= 测试（RWTH + 物理/模型分离） =========================
[x_input_test, fs_file] = audioread('bandpassed_200_700.wav');
x_in = resample(x_input_test(:,1), fs, fs_file);

[P_test, S_true, info_test] = get_pair(TEST_ROW);   % S_true = 真实物理次级路径

% ---- Phase A：用 Sc 在该测试任务上做一次辨识，得到 S_task_test ----
U_id_t = mk_noise(T_noise);
Y_id_t = filter(S_true,1,U_id_t);
cutSt  = randi([Ls_S, length(Y_id_t)]);
U_t    = U_id_t(cutSt-Ls_S+1 : cutSt);
Y_t    = Y_id_t(cutSt-Ls_S+1 : cutSt);

aS_test = MAML_Nstep_forget_S2(Ls_S);  % 你的类名如果是 _S 就用 _S
aS_test.Psi = Sc;                      % 用元初始化 Sc 起步
[~, S_task_test, ~] = aS_test.MAML_initial_S(U_t, Y_t, muS, lamdaS, epslonS);

% ---- 构造扰动 d（物理世界：P_test）----
Dis_1 = filter(P_test,1,x_in);

% ---- 用 FxLMS_phys 分离 S_true / S_model 进行测试 ----
Wc_zero = zeros(Len_N,1);
muw     = 5e-4;   % 与你原测试一致

% ① 模型=Sc（仅用元初始化S）
[e_zero_metaS, ~] = FxLMS_phys(Len_N, Wc_zero, Dis_1, x_in, S_true, Sc,           muw);
[e_maml_metaS, ~] = FxLMS_phys(Len_N, Wc,      Dis_1, x_in, S_true, Sc,           muw);

% ② 模型=S_task_test（Phase A 适配后的S）
[e_zero_adaptS, ~] = FxLMS_phys(Len_N, Wc_zero, Dis_1, x_in, S_true, S_task_test, muw);
[e_maml_adaptS, ~] = FxLMS_phys(Len_N, Wc,      Dis_1, x_in, S_true, S_task_test, muw);

% ③ 参考上界：模型=S_true（oracle）
[e_maml_oracleS, ~] = FxLMS_phys(Len_N, Wc,      Dis_1, x_in, S_true, S_true,      muw);

% ---- 可视化 ----
t = (0:length(Dis_1)-1)/fs;
figure;
plot(t, Dis_1, ...
     t, e_zero_metaS,  '--', ...
     t, e_maml_metaS,  '-',  ...
     t, e_zero_adaptS, '-.', ...
     t, e_maml_adaptS, '-',  ...
     t, e_maml_oracleS, ':', 'LineWidth',1.2);
grid on; xlabel('Time (s)'); ylabel('Error signal');
title(sprintf('Dual-meta test (row %d, sec=%d, prim=%d)', ...
      info_test.row, info_test.sec_col, info_test.prim_col));
legend({'ANC off',...
        'Zero-init W | S=Sc',...
        'Meta-init W | S=Sc',...
        'Zero-init W | S=S_{task}',...
        'Meta-init W | S=S_{task}',...
        'Meta-init W | S=S_{true} (oracle)'}, 'Location','best');

% 0) 先计算所有序列的最短长度，确保可裁剪
Llist = [numel(Dis_1), numel(e_zero_metaS), numel(e_maml_metaS), ...
         numel(e_zero_adaptS), numel(e_maml_adaptS), numel(e_maml_oracleS)];
Lmin  = min(Llist);

% 1) 去瞬态：头尾各裁去 M 点（取控制器/路径长度最大者，再与 Lmin 自适应收缩）
M0 = max([Len_N, numel(P_test), numel(S_true), numel(Sc), numel(S_task_test)]);
M  = min(M0, floor((Lmin-1)/2)-1);             % 防止越界
crop = @(x) x(M+1 : numel(x)-M);                % 裁后 t=0 已进入稳态区

d0            = crop(Dis_1(:));
e_z_metaS     = crop(e_zero_metaS(:));
e_w_metaS     = crop(e_maml_metaS(:));
e_z_adaptS    = crop(e_zero_adaptS(:));
e_w_adaptS    = crop(e_maml_adaptS(:));
e_w_oracleS   = crop(e_maml_oracleS(:));

% 2) 对称滑动均值（不因果），端点 'shrink'
win = 4096;                                     % 如需按秒设窗：win = round(0.2*fs);
ms_off       = movmean(d0.^2,          win, 'Endpoints','shrink');
ms_z_metaS   = movmean(e_z_metaS.^2,   win, 'Endpoints','shrink');
ms_w_metaS   = movmean(e_w_metaS.^2,   win, 'Endpoints','shrink');
ms_z_adaptS  = movmean(e_z_adaptS.^2,  win, 'Endpoints','shrink');
ms_w_adaptS  = movmean(e_w_adaptS.^2,  win, 'Endpoints','shrink');
ms_w_oracleS = movmean(e_w_oracleS.^2, win, 'Endpoints','shrink');

% 3) 时间轴（裁剪后从0开始）
tt = (0:numel(d0)-1)/fs;

% 4) 绘图
figure;
plot(tt,10*log10(ms_off       + 1e-10),'k' , 'LineWidth',1.8); hold on;
plot(tt,10*log10(ms_z_metaS   + 1e-10),'--', 'LineWidth',1.6);
plot(tt,10*log10(ms_w_metaS   + 1e-10),'-' , 'LineWidth',1.6);
plot(tt,10*log10(ms_z_adaptS  + 1e-10),':',  'LineWidth',1.6);
plot(tt,10*log10(ms_w_adaptS  + 1e-10),'-' , 'LineWidth',1.8);
plot(tt,10*log10(ms_w_oracleS + 1e-10),'-.', 'LineWidth',1.8);
grid on; xlabel('Time (s)'); ylabel('MSE (dB)');
title('Sliding-MSE (sym. window, shrink endpoints)');
legend({'ANC off', ...
        'Zero-init W | S=Sc', ...
        'Meta-init W | S=Sc', ...
        'Zero-init W | S=S_{task}', ...
        'Meta-init W | S=S_{task}', ...
        'Meta-init W | S=S_{true} (oracle)'}, 'Location','best');


%% ===== 对比次级路径：Sc vs 训练路径/测试路径/适配S_task（时域+频域+相位） =====
% 取两条训练任务对应的 S
[~, s_train1] = get_pair(TRAIN_ROWS(1));
[~, s_train2] = get_pair(TRAIN_ROWS(2));

% 统一长度到 Ls_S，并单位范数（仅比较“形状”）
S1    = trim_to_L(s_train1,   Ls_S);  S1    = S1   / max(norm(S1 ,2), 1e-12);
S2    = trim_to_L(s_train2,   Ls_S);  S2    = S2   / max(norm(S2 ,2), 1e-12);
Sc_n  = trim_to_L(Sc,         Ls_S);  Sc_n  = Sc_n / max(norm(Sc_n,2), 1e-12);
St_n  = trim_to_L(S_true,     Ls_S);  St_n  = St_n / max(norm(St_n,2), 1e-12);   % 测试真实路径
Stask = trim_to_L(S_task_test,Ls_S);  Stask = Stask/ max(norm(Stask,2), 1e-12);  % 测试估计路径

% —— 时域（前 256 taps）——
Lplot = min(256, Ls_S);
n_ms  = (0:Lplot-1)/fs*1e3;   % ms
figure;
subplot(2,1,1);
plot(n_ms, Sc_n (1:Lplot),'LineWidth',1.8); hold on;
plot(n_ms, S1   (1:Lplot),'--','LineWidth',1.2);
plot(n_ms, S2   (1:Lplot),':','LineWidth',1.2);
plot(n_ms, St_n (1:Lplot),'-.' ,'LineWidth',1.2);
plot(n_ms, Stask(1:Lplot),'-'  ,'LineWidth',1.4);
grid on; xlabel('Time (ms)'); ylabel('Amplitude (norm.)');
title('Secondary path (time domain, normalized)');
legend({'S_c (learned)','S_{train,11}','S_{train,23}','S_{test true}','S_{task est}'}, 'Location','best');

% —— 频域幅度（dB）——
Nfft = 4096;
[Hc,  fHz] = freqz(Sc_n ,1,Nfft,fs);
[H1,  ~ ]  = freqz(S1   ,1,Nfft,fs);
[H2,  ~ ]  = freqz(S2   ,1,Nfft,fs);
[Ht,  ~ ]  = freqz(St_n ,1,Nfft,fs);    % 真实测试路径
[He,  ~ ]  = freqz(Stask,1,Nfft,fs);    % 估计（适配）路径
% subplot(2,1,2);
figure;
plot(fHz, 20*log10(abs(Hc)+1e-12)-10,'LineWidth',1.8); hold on;
plot(fHz, 20*log10(abs(H1)+1e-12),'--','LineWidth',1.2);
plot(fHz, 20*log10(abs(H2)+1e-12),':','LineWidth',1.2);
plot(fHz, 20*log10(abs(Ht)+1e-12),'-.' ,'LineWidth',1.2);
plot(fHz, 20*log10(abs(He)+1e-12),'-'  ,'LineWidth',1.4);
xlim([0, fs/2]); grid on;
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB, norm.)');
title('Secondary path magnitude response (normalized)');
legend({'S_c (learned)','S_{train,11}','S_{train,23}','S_{test true}','S_{task est}'}, 'Location','best');

%% —— 相位比较（去整体时延后的相位残差，相对 S_{test true}）——
band = (fHz >= 200) & (fHz <= 1000);      % 只看训练/测试带
w = 2*pi*fHz/fs;                          % rad/sample

phi_rel_c = unwrap(angle(Hc./Ht));
phi_rel_1 = unwrap(angle(H1./Ht));
phi_rel_2 = unwrap(angle(H2./Ht));
phi_rel_e = unwrap(angle(He./Ht));

fit_c = polyfit(w(band), phi_rel_c(band), 1);   tau_c_samp = -fit_c(1);
fit_1 = polyfit(w(band), phi_rel_1(band), 1);   tau_1_samp = -fit_1(1);
fit_2 = polyfit(w(band), phi_rel_2(band), 1);   tau_2_samp = -fit_2(1);
fit_e = polyfit(w(band), phi_rel_e(band), 1);   tau_e_samp = -fit_e(1);

phi_res_c = phi_rel_c - polyval(fit_c, w);
phi_res_1 = phi_rel_1 - polyval(fit_1, w);
phi_res_2 = phi_rel_2 - polyval(fit_2, w);
phi_res_e = phi_rel_e - polyval(fit_e, w);

figure;
plot(fHz, rad2deg(phi_res_c), 'LineWidth', 1.8); hold on;
plot(fHz, rad2deg(phi_res_1), '--', 'LineWidth', 1.2);
plot(fHz, rad2deg(phi_res_2), ':',  'LineWidth', 1.2);
plot(fHz, rad2deg(phi_res_e), '-' , 'LineWidth', 1.4);
xlim([0, fs/2]); grid on;
xlabel('Frequency (Hz)'); ylabel('Phase residual (deg)');
title('Secondary path phase (delay-neutral, relative to S_{test true})');
legend({'S_c vs S_{test}','S_{train,11} vs S_{test}','S_{train,23} vs S_{test}','S_{task} vs S_{test}'}, 'Location','best');

% —— 定量指标（200–1000 Hz）——
mag_rel = @(Ha,Hb) norm(abs(Ha(band))-abs(Hb(band))) / max(norm(abs(Hb(band))),1e-12);
coh     = @(Ha,Hb) (Ha(band)'*conj(Hb(band))) / sqrt( (Ha(band)'*conj(Ha(band))) * (Hb(band)'*conj(Hb(band))) );

mag_err_c = mag_rel(Hc,Ht);   ph_std_c = std(rad2deg(phi_res_c(band)));  rho_c = abs(coh(Hc,Ht));
mag_err_1 = mag_rel(H1,Ht);   ph_std_1 = std(rad2deg(phi_res_1(band)));  rho_1 = abs(coh(H1,Ht));
mag_err_2 = mag_rel(H2,Ht);   ph_std_2 = std(rad2deg(phi_res_2(band)));  rho_2 = abs(coh(H2,Ht));
mag_err_e = mag_rel(He,Ht);   ph_std_e = std(rad2deg(phi_res_e(band)));  rho_e = abs(coh(He,Ht));

fprintf(['[S vs S_{test true} | 200–1000 Hz]\n' ...
         '  S_c   : mag_err=%.2f%%, tau=%.2f samp (%.3f ms), phase_std=%.2f deg, |rho|=%.3f\n' ...
         '  S_11  : mag_err=%.2f%%, tau=%.2f samp (%.3f ms), phase_std=%.2f deg, |rho|=%.3f\n' ...
         '  S_23  : mag_err=%.2f%%, tau=%.2f samp (%.3f ms), phase_std=%.2f deg, |rho|=%.3f\n' ...
         '  S_task: mag_err=%.2f%%, tau=%.2f samp (%.3f ms), phase_std=%.2f deg, |rho|=%.3f\n'], ...
         100*mag_err_c, tau_c_samp, 1e3*tau_c_samp/fs, ph_std_c, rho_c, ...
         100*mag_err_1, tau_1_samp, 1e3*tau_1_samp/fs, ph_std_1, rho_1, ...
         100*mag_err_2, tau_2_samp, 1e3*tau_2_samp/fs, ph_std_2, rho_2, ...
         100*mag_err_e, tau_e_samp, 1e3*tau_e_samp/fs, ph_std_e, rho_e);

%% ============== 跨任务泛化对比：Task1 训练 → Task2 测试 ==============
TASK1_ROW = BASE_ROW;     % 训练任务（task1）
TASK2_ROW = TEST_ROW;     % 测试任务（task2）

% 取 task1 / task2 的路径
[P1, S1] = get_pair(TASK1_ROW);
[P2, S2] = get_pair(TASK2_ROW);

% ---------------- Task1：训练得到 W_task1 ----------------
N_train = 4 * length(x_in);                 % 训练时长（可调）
x_train = filter(bp,1, randn(N_train,1));   % 训练参考
Dis_tr1 = filter(P1,1,x_train);             % 物理扰动（task1）

muw_tr  = 0.03;                              % 训练期步长（沿用你原设置）
% 在 task1 上训练：S_true = S1, S_model = S1
[~, W_task1] = FxLMS_phys(Len_N, zeros(Len_N,1), Dis_tr1, x_train, S1, S1, muw_tr);

% ---------------- Task2：测试期（真实路径= S2） ----------------
Dis_t2 = filter(P2,1,x_in);                 % 物理扰动（task2）
muw_te = 5e-4;                              % 测试期步长

% 对比①：Transfer（W=W_task1, S_model=S1, S_true=S2）
[e_transfer, ~] = FxLMS_phys(Len_N, W_task1, Dis_t2, x_in, S2, S1, muw_te);

% 对比②：MAML（W=Wc, S_model=Sc, S_true=S2）
[e_maml_xfer, ~] = FxLMS_phys(Len_N, Wc,      Dis_t2, x_in, S2, Sc, muw_te);

% 可视化（时域误差）
t = (0:length(Dis_t2)-1)/fs;
figure;
plot(t, Dis_t2, 'k', ...
     t, e_transfer,  '--', ...
     t, e_maml_xfer, '-', 'LineWidth',1.2);
grid on; xlabel('Time (s)'); ylabel('Error signal');
title(sprintf('Cross-task: train row %d → test row %d', TASK1_ROW, TASK2_ROW));
legend({'ANC off', 'Transfer: W_{task1}, S_{model}=S_1', 'MAML: W_c, S_{model}=S_c'}, 'Location','best');

%% ---------- MSE（去瞬态 + 对称窗） ----------
Llist = [numel(Dis_t2), numel(e_transfer), numel(e_maml_xfer)];
Lmin  = min(Llist);

M0 = max([Len_N, numel(P2), numel(S2), numel(S1), numel(Sc)]);
M  = min(M0, floor((Lmin-1)/2)-1);                % 防越界
crop = @(x) x(M+1 : numel(x)-M);

d0        = crop(Dis_t2(:));
e_tr_c    = crop(e_transfer(:));
e_ma_c    = crop(e_maml_xfer(:));

win = 4096;
ms_off  = movmean(d0   .^2, win, 'Endpoints','shrink');
ms_tr   = movmean(e_tr_c.^2, win, 'Endpoints','shrink');
ms_maml = movmean(e_ma_c.^2, win, 'Endpoints','shrink');

tt = (0:numel(d0)-1)/fs;
figure;
plot(tt,10*log10(ms_off +1e-10),'k','LineWidth',1.8); hold on;
plot(tt,10*log10(ms_tr  +1e-10),'--','LineWidth',1.6);
plot(tt,10*log10(ms_maml+1e-10),'-' ,'LineWidth',1.8);
grid on; xlabel('Time (s)'); ylabel('MSE (dB)');
title('Cross-task Sliding-MSE (sym. window, shrink endpoints)');
legend({'ANC off','Transfer: W_{task1}, S_{model}=S_1','MAML: W_c, S_{model}=S_c'}, 'Location','best');
%% ====================（可选）简单 MSE 对比图（from t=0, movmean）====================
% 去瞬态：裁掉控制滤波器/主路径/次级路径中的最长长度（只裁开头，保留全长尾部）
% M    = max([Len_N, numel(P_test), numel(S_test)]);
% crop = @(x) x(M+1:end);
% 
% d0       = crop(Dis_1(:));
% e_zero   = crop(Er_zero(:));
% e_single = crop(Er_single(:));
% e_maml   = crop(Er_maml(:));
% 
% % 统一到最短长度，保证时间轴一致
% L = min([numel(d0), numel(e_zero), numel(e_single), numel(e_maml)]);
% d0       = d0(1:L);
% e_zero   = e_zero(1:L);
% e_single = e_single(1:L);
% e_maml   = e_maml(1:L);
% 
% % ===== 因果滑动均值 MSE：第 n 点只平均 [n-win+1, n]（起点自动用短窗）=====
% win   = 4096;                                 % 可调短些看更早期动态
% count = filter(ones(win,1), 1, ones(L,1));    % 每点的实际窗长（起点<win）
% pow   = @(x) filter(ones(win,1), 1, x(:).^2) ./ count;   % 因果均值功率
% toDb  = @(x) 10*log10(x + 1e-10);
% 
% ms_off    = toDb(pow(d0));
% ms_zero   = toDb(pow(e_zero));
% ms_single = toDb(pow(e_single));
% ms_maml   = toDb(pow(e_maml));
% 
% tt = (0:L-1)/fs;
% 
% figure;
% plot(tt, ms_off,    'k',  'LineWidth',1.8); hold on;
% plot(tt, ms_zero,   '--', 'LineWidth',1.6);
% plot(tt, ms_single, '-.', 'LineWidth',1.6);
% plot(tt, ms_maml,   '-',  'LineWidth',1.8);
% xlabel('Time (s)'); ylabel('MSE (dB)'); grid on;
% title('Sliding-MSE (causal, starts at t=0)');
% legend({'ANC off','Zero init','Single-task init','MAML init'}, 'Location','best');

%% ======================= 本脚本用到的小工具 =======================
function [P,S,info] = local_pick_pair(Hpri,Hsec,meta,row,fs,Fs0)
    % 由 meta 的“行号”拿到成对列号，再取列并重采样
    sc = meta.sec_col(row);    % 对应次级列
    pc = meta.prim_col(row);   % 对应初级列
    P  = Hpri(:,pc);  S = Hsec(:,sc);
    if fs ~= Fs0
        P = resample(P, fs, Fs0);
        S = resample(S, fs, Fs0);
    end
    P = P(:); S = S(:);
    info.row = row;
    info.sec_col = sc;
    info.prim_col = pc;
end

function hL = trim_to_L(h, L)
    h = h(:);
    if numel(h) >= L
        hL = h(1:L);
    else
        hL = [h; zeros(L-numel(h),1)];
    end
end



