function theta = run_ls_filter(x_s, xd_s, u_s, Ts, lambda)
% Filter-based LS helper extracted from part2.m
% Inputs: x_s (position samples), xd_s (velocity samples), u_s (input samples), Ts (sampling), lambda (filter pole)
% Returns: theta = [m; c; k]

    % Bilinear coefficients for H1(z) = 1/(s+lambda)
    b1 = Ts * [1,  1];
    a1 = [2 + lambda*Ts,  lambda*Ts - 2];
    % Bilinear coefficients for Hs(z) = s/(s+lambda)
    bs = 2  * [1, -1];
    % (same denominator a1)

    % Filtered signals
    z     = filter(b1, a1, u_s);    % (1/Lambda) * u
    zeta1 = filter(bs, a1, xd_s);  % (s/Lambda) * xd  [replaces xdd/Lambda]
    zeta2 = filter(b1, a1, xd_s);  % (1/Lambda) * xd
    zeta3 = filter(b1, a1, x_s);   % (1/Lambda) * x

    % Skip filter settling transient (~5 time constants = 5/lambda seconds)
    skip = ceil(5 / (lambda * Ts));
    ii   = (skip + 1) : length(z);
    if length(ii) < 10
        theta = nan(3, 1);
        return;
    end

    Z     = [zeta1(ii)', zeta2(ii)', zeta3(ii)'];
    z_vec = z(ii)';
    theta = (Z' * Z) \ (Z' * z_vec);
end