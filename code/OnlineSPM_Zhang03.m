% Copyright (c) 2026 Ziyi Yang (ziyi016@e.ntu.edu.sg).
% Released for the ICASSP 2026 MAML co-initialization experiments.
%
function [S_hat, out] = OnlineSPM_Zhang03( ...
    x_ref, d_mic, S_true, init, ...
    Lw, Ls, Lh, mu_w, mu_s, mu_h, varargin)
% Robust online secondary-path modeling (Zhang, Lan, Ser, 2003)
% Implements the OSPM-FxLMS flow used in the paper:
% e = d - S*u + S*vm, e1^pi = e - v_hat, the secondary-path
% estimate is updated with vms*es, and the control-filter update uses the
% current estimated secondary path for filtered-x.
%
% Required:
%   x_ref(Nx1)  : reference signal.
%   d_mic(Nx1)  : disturbance at the error microphone before control/probe.
%   S_true(Lsx1): true secondary path.
%   init        : struct with init.w0(Lw), init.s0(Ls), init.h0(Lh).
%   Lw,Ls,Lh    : filter lengths.
%   mu_w,mu_s,mu_h : step sizes.
%
% Name-Value:
%   'Tw','Ts','Th'            : norm constraints, default 10.
%   'alpha_power'             : power EMA factor, default 0.99.
%   'c_aux'                   : auxiliary-noise power constant, default 1.0.
%   'delay'                   : integer reference-to-error delay, default 0.
%   'SNR_ref_dB','SNR_err_dB' : optional measured AWGN on x_ref / d_mic.

% ---------- Parse options ----------
opt = struct('Tw',10,'Ts',10,'Th',10, ...
             'alpha_power',0.99,'c_aux',1.0,'delay',0, ...
             'SNR_ref_dB',[],'SNR_err_dB',[]);
for k = 1:2:numel(varargin)
    opt.(varargin{k}) = varargin{k+1};
end

% ---------- Preprocess and initialize ----------
x_ref = x_ref(:); d_mic = d_mic(:); S_true = S_true(:);
N = min(numel(x_ref), numel(d_mic));
x_ref = x_ref(1:N); d_mic = d_mic(1:N);

if ~isempty(opt.SNR_ref_dB), x_ref = add_awgn_measured_local(x_ref, opt.SNR_ref_dB); end
if ~isempty(opt.SNR_err_dB), d_mic = add_awgn_measured_local(d_mic, opt.SNR_err_dB); end

w  = zeros(Lw,1); if isfield(init,'w0') && ~isempty(init.w0), w  = init.w0(:); end
sp = zeros(Ls,1); if isfield(init,'s0') && ~isempty(init.s0), sp = init.s0(:); end
h  = zeros(Lh,1); if isfield(init,'h0') && ~isempty(init.h0), h  = init.h0(:); end

% Shift-register buffers.
xw   = zeros(Lw,1);                % Reference taps for w.
xh   = zeros(Lh,1);                % Reference taps for h.
xs   = zeros(Ls+Lw,1);             % Buffer for filtered-x synthesis.
f_w1 = zeros(Ls,1);                % Control-signal history.
vms  = zeros(Ls,1);                % Auxiliary-noise history.

Px = 1; Pe1_pie = 0;

% Logs.
logS = zeros(Ls,N); logW = zeros(Lw,N); logH = zeros(Lh,N);
e = zeros(N,1); e_no_aux = zeros(N,1); e1_pie = zeros(N,1); e1h = zeros(N,1);
aux_at_mic = zeros(N,1);

% ---------- Main loop ----------
for n = 1:N
    % taps
    xw = [x_ref(n); xw(1:end-1)];
    if opt.delay>0
        if n>opt.delay, xh = [x_ref(n-opt.delay); xh(1:end-1)];
        else,           xh = [0;                 xh(1:end-1)];
        end
    else
        xh = [x_ref(n); xh(1:end-1)];
    end
    xs = [xs(2:end); x_ref(n)];

    % Control drive u(n).
    u = w.'*xw;

    % Auxiliary-noise power schedule.
    v  = randn();
    if Pe1_pie < Px
        vm = opt.c_aux * v * sqrt(Pe1_pie + 1e-12);
    else
        vm = opt.c_aux * v * sqrt(Px + 1e-12);
    end

    f_w1 = [u;  f_w1(1:end-1)];
    vms  = [vm; vms(1:end-1)];

    % Physical response at the microphone.
    y_u   = f_w1.' * S_true;           % S * u
    y_vm  = vms.'  * S_true;           % S * vm
    e(n)  = d_mic(n) - y_u + y_vm;     % e = d - S*u + S*vm
    e_no_aux(n) = d_mic(n) - y_u;

    % Auxiliary estimate using the current secondary-path model.
    v_hat = vms.' * sp;

    % Residuals for w/h updates.
    e1_pie(n) = e(n) - v_hat;
    z_hat     = xh.' * h;
    e1h(n)    = e1_pie(n) - z_hat;

    % Residual for the secondary-path update.
    g  = e(n) - z_hat;
    es = g - v_hat;

    % 1) Estimated secondary path.
    sp_new = sp + mu_s * (vms * es);
    if norm(sp_new) <= opt.Ts
        sp = sp_new;
    else
        sp = opt.Ts * sp_new / (norm(sp_new)+1e-12);
    end

    % 2) Control filter using current estimated secondary path.
    temp_fx = filter(sp,1,xs);
    f_s11_  = temp_fx(end:-1:end-Lw+1);
    w_new   = w + mu_w * f_s11_ * e1_pie(n);
    if norm(w_new) <= opt.Tw
        w = w_new;
    else
        w = opt.Tw * w_new / (norm(w_new)+1e-12);
    end

    % 3) Reference-to-error canceller.
    h_new = h + mu_h * xh * e1h(n);
    if norm(h_new) <= opt.Th
        h = h_new;
    else
        h = opt.Th * h_new / (norm(h_new)+1e-12);
    end

    % Power EMAs.
    Px      = opt.alpha_power * Px      + (1-opt.alpha_power) * (x_ref(n)^2);
    Pe1_pie = opt.alpha_power * Pe1_pie + (1-opt.alpha_power) * (e1_pie(n)^2);

    % Logs.
    logS(:,n)=sp; logW(:,n)=w; logH(:,n)=h;
    aux_at_mic(n) = y_vm;
end

% ---------- Outputs ----------
S_hat = sp;
out.w_end = w; out.h_end = h;
out.S_trace = logS; out.W_trace = logW; out.H_trace = logH;
out.e = e; out.e_no_aux = e_no_aux; out.e1_pie = e1_pie; out.e1h = e1h;
out.aux_at_mic = aux_at_mic; out.Px = Px; out.Pe1_pie = Pe1_pie;
end

function y = add_awgn_measured_local(x, snr_db)
    p_signal = mean(x(:).^2);
    p_noise = p_signal / (10^(snr_db/10));
    y = x + sqrt(p_noise) * randn(size(x));
end
