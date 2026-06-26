% Copyright (c) 2026 Ziyi Yang (ziyi016@e.ntu.edu.sg).
% Extended for the ICASSP 2026 MAML co-initialization experiments.
% Portions of the FxLMS code lineage trace to Shi Dongyuan's Meta repository:
% https://github.com/ShiDongyuan/Meta
%
function [e, w] = FxLMS_phys(L, w0, d, x, S_true, S_model, mu)
% L: control-filter length.
% w0: initial control filter.
% d: disturbance at the error microphone.
% x: raw reference signal.
% S_true: physical secondary path used to synthesize the microphone error.
% S_model: estimated secondary path used to build the filtered-x regressor.
    w = w0(:); x = x(:); d = d(:);
    N = length(d);
    e = zeros(N,1);

    % Precompute the model-based filtered reference.
    xprime = filter(S_model,1,x);

    % Convolve the loudspeaker drive with the true secondary path.
    Ls = length(S_true); buf = zeros(Ls,1);

    for n = L:N
        x_tap  = x(n:-1:n-L+1);           % Raw reference taps.
        xp_tap = xprime(n:-1:n-L+1);      % Filtered-x taps for the gradient.

        u = w.'*x_tap;                    % Loudspeaker drive.
        buf = [u; buf(1:end-1)];
        y = S_true.'*buf;                 % Physical anti-noise output.

        e(n) = d(n) - y;

        den = xp_tap.'*xp_tap + 1e-8;     % NLMS
        w   = w + mu * e(n) * (xp_tap / den);
    end
end
