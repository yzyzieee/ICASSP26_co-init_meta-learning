% OnlineSPM_ThreePhase_RWTH_switchS.m
% 三阶段连续收敛：Phase1 用 (P1,S1)；Phase2 切到 (P2,S2)；Phase3 再切到 (P3,S3)
% 传统：后续相位沿用上一相位的 w,s,h
% 提出：每次切换把 w,s 重置为 (Wc,Sc)
% 在线建模内核完全照“程序A”(Zhang03) 的公式/流程

clear; clc;

%% ============ 基本配置 ============
fs = 16000;
PAIRS_FILE = 'ANC_PriSec_pairs.mat';
TASK1_ROW  = 10;      % Phase1 的 (P1,S1)
TASK2_ROW  = 25;      % Phase2 的 (P2,S2)
TASK3_ROW  = 38;      % === Phase3 的 (P3,S3)（可按需要修改为有效行号） ===

band_hz = [200 1000]; fir_ord = 512;

% 程序A一致的长度/参数
N=512; L_h=512; L_s=512; delay=0;
SNR=30;
T1=60; T2=60; T3=60;      % === 三个相位均为 60 s ===
it1=T1*fs; it2=T2*fs; it3=T3*fs;

alpha_power = 0.99;
c_aux = 1.0;
mu_w = 1e-2;  mu_s = 1e-3;  mu_h = 1e-3;
Tw=10; Ts=10; Th=10;

% 读取 MAML 初始化
% Sws = load('Wc_Sc_MAML_RWTH.mat','Wc','Sc');
Sws = load('Wc_Sc_MAML_RWTH_Min_3.mat','Wc','Sc');
Wc0 = pad_or_trim(Sws.Wc(:), N);
Sc0 = pad_or_trim(Sws.Sc(:), L_s);

%% ============ 载入三组 (P,S) ============
S0   = load(PAIRS_FILE,'out_pairs');
Fs0  = S0.out_pairs.Fs;
Hpri = S0.out_pairs.H_primary;
Hsec = S0.out_pairs.H_secondary;
meta = S0.out_pairs.meta_pairs;

% row -> (P,S)
[P1,S1] = pick_pair(meta,Hpri,Hsec,TASK1_ROW,fs,Fs0);
[P2,S2] = pick_pair(meta,Hpri,Hsec,TASK2_ROW,fs,Fs0);
[P3,S3] = pick_pair(meta,Hpri,Hsec,TASK3_ROW,fs,Fs0);   % === Phase3 ===
p1 = pad_or_trim(P1(:),N);  s1 = pad_or_trim(S1(:),L_s);
p2 = pad_or_trim(P2(:),N);  s2 = pad_or_trim(S2(:),L_s);
p3 = pad_or_trim(P3(:),N);  s3 = pad_or_trim(S3(:),L_s);

%% ============ 构造整段参考 & 扰动 ============
bp  = fir1(fir_ord, band_hz/(fs/2));
ref = filter(bp,1, randn(it1+it2+it3,1));      % 带限参考（总长度三段相加）
ref_noisy = awgn(ref, SNR, 'measured');        % 程序A：参考加测量噪声

% 三段使用各自的 Pk 生成 d
d1 = filter(p1,1, ref(1:it1));
d2 = filter(p2,1, ref(it1+1:it1+it2));
d3 = filter(p3,1, ref(it1+it2+1:end));
d1n = awgn(d1, SNR, 'measured');
d2n = awgn(d2, SNR, 'measured');
d3n = awgn(d3, SNR, 'measured');

ref1 = ref_noisy(1:it1);
ref2 = ref_noisy(it1+1:it1+it2);
ref3 = ref_noisy(it1+it2+1:end);               % === Phase3 ===

% —— Phase 2 抬高 +3 dB（功率），即振幅乘以 10^(3/20)
boost_dB_P2 = 2.5; g2 = 10^(boost_dB_P2/20);
ref2 = ref2 * g2;  d2n = d2n * g2;

% （如需 Phase 3 也抬升，可取消注释）
boost_dB_P3 = 3.5; g3 = 10^(boost_dB_P3/20);
ref3 = ref3 * g3;  d3n = d3n * g3;

%% ============ Phase1：Zero-init vs MAML-init（真路径=S1） ============
rng(1234);
[eZ1, stZ1] = run_A_phase(ref1, d1n, s1, ...
    zeros(N,1), zeros(L_s,1), zeros(L_h,1), ...
    mu_w,mu_s,mu_h, Tw,Ts,Th, alpha_power,c_aux, delay);

rng(1234);
[eM1, stM1] = run_A_phase(ref1, d1n, s1, ...
    Wc0,       Sc0,        zeros(L_h,1), ...
    mu_w,mu_s,mu_h, Tw,Ts,Th, alpha_power,c_aux, delay);

%% ============ Phase2：切换真实路径=S2 & (P2) ============
% 传统：沿用 Phase1 的 (w,s,h)
% 提出：把 (w,s) 重置为 (Wc0,Sc0)
rng(5678);
[eZ2, stZ2] = run_A_phase(ref2, d2n, s2, ...
    stZ1.w,    stZ1.s,     stZ1.h, ...
    mu_w,mu_s,mu_h, Tw,Ts,Th, alpha_power,c_aux, delay);

rng(5678);
[eM2, stM2] = run_A_phase(ref2, d2n, s2, ...
    Wc0,       Sc0,        zeros(L_h,1), ...
    mu_w,mu_s,mu_h, Tw,Ts,Th, alpha_power,c_aux, delay);

%% ============ Phase3：再切换真实路径=S3 & (P3) ============
% 传统：沿用 Phase2 的 (w,s,h)
% 提出：再次在切换处把 (w,s) 重置为 (Wc0,Sc0)
rng(9012);
[eZ3, stZ3] = run_A_phase(ref3, d3n, s3, ...
    stZ2.w,    stZ2.s,     stZ2.h, ...
    mu_w,mu_s,mu_h, Tw,Ts,Th, alpha_power,c_aux, delay);

rng(9012);
[eM3, stM3] = run_A_phase(ref3, d3n, s3, ...
    Wc0,       Sc0,        zeros(L_h,1), ...
    mu_w,mu_s,mu_h, Tw,Ts,Th, alpha_power,c_aux, delay);

%% ============ 拼接三阶段 ============
e_zero = [eZ1; eZ2; eZ3];
e_maml = [eM1; eM2; eM3];
d_all  = [d1n; d2n; d3n];

% === 三阶段拼接后的辅助噪声（@mic） ===
auxZ_all = [stZ1.aux; stZ2.aux; stZ3.aux];
auxM_all = [stM1.aux; stM2.aux; stM3.aux];

%% ============ Sliding-MSE（三阶段连续曲线） ============
win = 4096;
ms_off = movmean(d_all.^2 , win, 'Endpoints','shrink');
ms_Z   = movmean(e_zero.^2, win, 'Endpoints','shrink');
ms_M   = movmean(e_maml.^2, win, 'Endpoints','shrink');

tt    = (0:numel(d_all)-1)/fs;
t_bd1 = T1;
t_bd2 = T1 + T2;   % === 第二个切换点 ===

figure; hold on; grid on;
plot(tt,10*log10(ms_off+1e-10),'k','LineWidth',1.6);
plot(tt,10*log10(ms_Z  +1e-10),'--','LineWidth',1.6);
plot(tt,10*log10(ms_M  +1e-10),'-','LineWidth',1.8);
xline(5,'k:','Phase 1 (S_1, P_1)', ...
    'LabelOrientation','horizontal','LabelVerticalAlignment','bottom');
xline(t_bd1,'k:','Phase 2 (switch to S_2, P_2)', ...
    'LabelOrientation','horizontal','LabelVerticalAlignment','bottom');
xline(t_bd2,'k:','Phase 3 (switch to S_3, P_3)', ...
    'LabelOrientation','horizontal','LabelVerticalAlignment','bottom');
xlabel('Time (s)'); ylabel('MSE (dB)');
title('Three-phase online modeling: real S switches at t = 60 s and 120 s');
legend({'ANC off','Zero-init','MAML-init (reset at phase boundaries)'}, 'Location','best');

%% ============ 误差 + 辅助噪声功率（同窗） ============
ms_aux_Z = movmean(auxZ_all.^2, win, 'Endpoints','shrink');
ms_aux_M = movmean(auxM_all.^2, win, 'Endpoints','shrink');

figure; hold on; grid on;
plot(tt, 10*log10(ms_off+1e-12), 'k' , 'LineWidth',1.6);         % ANC off
plot(tt, 10*log10(ms_Z  +1e-12), '--', 'LineWidth',1.6);         % Zero-init 误差
plot(tt, 10*log10(ms_M  +1e-12), '-' , 'LineWidth',1.8);         % MAML-init 误差
plot(tt, 10*log10(ms_aux_Z+1e-12), ':', 'LineWidth',1.6, 'Color',[0.65 0.20 0.80]); % Zero-init 辅助
plot(tt, 10*log10(ms_aux_M+1e-12), ':', 'LineWidth',1.6, 'Color',[0.20 0.35 0.85]); % MAML-init 辅助
xline(5,'k:','Phase 1 (S_1, P_1)', ...
    'LabelOrientation','horizontal','LabelVerticalAlignment','bottom');
xline(t_bd1,'k:','Phase 2 (switch to S_2, P_2)', ...
      'LabelOrientation','horizontal','LabelVerticalAlignment','bottom');
xline(t_bd2,'k:','Phase 3 (switch to S_3, P_3)', ...
      'LabelOrientation','horizontal','LabelVerticalAlignment','bottom');
xlabel('Time (s)'); ylabel('Mean Squared Error (dB)');
title('Online modeling FxLMS with auxiliary-noise power (switches at 60 s & 120 s)');
legend({'ANC off', 'Error (Original)', 'Error (Proposed MAML)', ...
        'Aux noise (Original)', 'Aux noise (Proposed MAML)'}, ...
        'Location','best');


%% ===== 子程序 =====
function [e, st] = run_A_phase(ref_noisy, d_noisy, S_true, ...
                               w0, s0, h0, ...
                               mu_w,mu_s,mu_h, Tw,Ts,Th, alpha_power,c_aux, delay)
    N=numel(w0); L_s=numel(s0); L_h=numel(h0); it=numel(ref_noisy);
    e=zeros(it,1); xw=zeros(N,1); xh=zeros(L_h,1); xs=zeros(L_s+N,1);
    f_w1=zeros(L_s,1); vms=zeros(L_s,1);
    w=w0; sp_=s0; h=h0; Px=1; Pe1_pie=0;

    % ==== 新增：记录辅助噪声到麦克风 ====
    aux_mic = zeros(it,1);

    for n=1:it
        xw=[ref_noisy(n);xw(1:end-1)];
        if n<delay+1, xh=[0;xh(1:end-1)]; else, xh=[ref_noisy(n-delay);xh(1:end-1)]; end
        xs=[xs(2:end);ref_noisy(n)];
        f_w1=[w.'*xw; f_w1(1:end-1)];

        v=randn();
        if Pe1_pie<(Px/1), vm=c_aux*v*sqrt(Pe1_pie+1e-12);
        else,              vm=c_aux*v*sqrt(Px+1e-12);
        end
        vms=[vm; vms(1:end-1)];

        y_u = f_w1.'*S_true;
        v_m = vms.' *S_true;
        e(n)= d_noisy(n) - y_u + v_m;

        % ==== 新增：保存到数组 ====
        aux_mic(n) = v_m;

        v_hat  = vms.'*sp_;
        e1_pie = e(n) - v_hat;
        z_hat  = xh.'*h;
        e1h    = e1_pie - z_hat;
        g      = e(n) - z_hat;
        es     = g - v_hat;

        sp_new = sp_ + mu_s*(vms*es);
        if norm(sp_new)<=Ts, sp_=sp_new; else, sp_=Ts*sp_new/(norm(sp_new)+1e-12); end

        tmp = filter(sp_,1,xs);
        f_s11_ = tmp(end:-1:end-N+1);

        w_new = w + mu_w*(f_s11_*e1_pie);
        if norm(w_new)<=Tw, w=w_new; else, w=Tw*w_new/(norm(w_new)+1e-12); end

        h_new = h + mu_h*(xh*e1h);
        if norm(h_new)<=Th, h=h_new; else, h=Th*h_new/(norm(h_new)+1e-12); end

        Px      = alpha_power*Px      + (1-alpha_power)*(ref_noisy(n)^2);
        Pe1_pie = alpha_power*Pe1_pie + (1-alpha_power)*(e1_pie^2);
    end

    st.w=w; st.s=sp_; st.h=h;
    % ==== 新增：返回辅助噪声轨迹 ====
    st.aux = aux_mic;
end

function [P,S] = pick_pair(meta,Hpri,Hsec,row,fs,Fs0)
    pc=meta.prim_col(row); sc=meta.sec_col(row);
    P=Hpri(:,pc); S=Hsec(:,sc);
    if fs~=Fs0, P=resample(P,fs,Fs0); S=resample(S,fs,Fs0); end
end

function xL = pad_or_trim(x,L)
    x=x(:); if numel(x)>=L, xL=x(1:L); else, xL=[x;zeros(L-numel(x),1)]; end
end
