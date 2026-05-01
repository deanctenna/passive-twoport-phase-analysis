function [lambda_star, x_star, history] = solve_min_normalized_real(A, opts)
%SOLVE_MIN_NORMALIZED_REAL Minimize Re(x'*A*x) / (||x||_2 ||A*x||_2) for 2x2 A.
%
% Usage:
%   [lambda_star, x_star, history] = solve_min_normalized_real(A)
%   [lambda_star, x_star, history] = solve_min_normalized_real(A, opts)
%
% Optional fields in opts:
%   tol         - outer stopping tolerance (default: 1e-8)
%   max_outer   - maximum outer Dinkelbach iterations (default: 50)
%   n_theta     - number of theta grid points (default: 121)
%   n_phi       - number of phi grid points (default: 241)
%   n_starts    - number of best grid points used for local refinement (default: 8)
%   ax_tol      - threshold for treating ||A*x||_2 as zero (default: 1e-12)
%   display     - true to print progress (default: true)
%
% history columns:
%   [iteration, lambda_k, residual_k, theta_k, phi_k]

if nargin < 1 || isempty(A)
    A = input('Please enter a 2x2 complex matrix A, for example [1 2+1i; -3i 4]: ');
end

if nargin < 2
    opts = struct();
end

if ~isequal(size(A), [2, 2])
    error('A must be a 2x2 complex matrix.');
end

tol = get_option(opts, 'tol', 1e-8);
max_outer = get_option(opts, 'max_outer', 50);
n_theta = get_option(opts, 'n_theta', 121);
n_phi = get_option(opts, 'n_phi', 241);
n_starts = get_option(opts, 'n_starts', 8);
ax_tol = get_option(opts, 'ax_tol', 1e-12);
display_progress = get_option(opts, 'display', true);

lambda_k = 0;
history = nan(max_outer, 5);
x_star = [];

for k = 1:max_outer
    [theta_k, phi_k] = solve_inner_problem(A, lambda_k, n_theta, n_phi, n_starts, ax_tol);
    x_k = unit_vector(theta_k, phi_k);
    Ax_k = A * x_k;
    gx_k = norm(Ax_k, 2);

    if gx_k <= ax_tol
        error('Encountered ||A*x||_2 approximately zero during iteration %d.', k);
    end

    fx_k = real(x_k' * A * x_k);
    residual_k = fx_k - lambda_k * gx_k;

    history(k, :) = [k - 1, lambda_k, residual_k, theta_k, phi_k];
    x_star = x_k;

    if display_progress
        fprintf('iter = %2d, lambda = %.12f, residual = %.3e\n', ...
            k - 1, lambda_k, residual_k);
    end

    if abs(residual_k) < tol
        lambda_star = lambda_k;
        history = history(1:k, :);
        return;
    end

    lambda_k = fx_k / gx_k;
end

lambda_star = lambda_k;
warning('Maximum outer iterations reached before the residual met the tolerance.');

end

function [theta_best, phi_best] = solve_inner_problem(A, lambda_k, n_theta, n_phi, n_starts, ax_tol)
theta_grid = linspace(0, pi / 2, n_theta);
phi_grid = linspace(0, 2 * pi, n_phi + 1);
phi_grid(end) = [];

values = inf(n_theta, n_phi);

for i = 1:n_theta
    theta = theta_grid(i);
    for j = 1:n_phi
        phi = phi_grid(j);
        values(i, j) = inner_objective(theta, phi, A, lambda_k, ax_tol);
    end
end

values_vector = values(:);
[~, order] = sort(values_vector, 'ascend');
num_candidates = min(n_starts, numel(order));

theta_best = theta_grid(1);
phi_best = phi_grid(1);
best_value = inf;

options = optimset( ...
    'Display', 'off', ...
    'MaxIter', 500, ...
    'MaxFunEvals', 1500, ...
    'TolX', 1e-10, ...
    'TolFun', 1e-10);

for idx = 1:num_candidates
    linear_index = order(idx);
     
    [i, j] = ind2sub(size(values), linear_index);

    theta0 = theta_grid(i);
    phi0 = phi_grid(j);
    uv0 = [theta_to_u(theta0); phi0];

    objective_uv = @(uv) inner_objective_uv(uv, A, lambda_k, ax_tol);
    [uv_opt, ~] = fminsearch(objective_uv, uv0, options);
    [theta_opt, phi_opt] = unpack_uv(uv_opt);
    local_value = inner_objective(theta_opt, phi_opt, A, lambda_k, ax_tol);

    if local_value < best_value
        best_value = local_value;
        theta_best = theta_opt;
        phi_best = phi_opt;
    end
end

if ~isfinite(best_value)
    error('The inner problem did not find a feasible point with ||A*x||_2 > 0.');
end

end

function value = inner_objective_uv(uv, A, lambda_k, ax_tol)
[theta, phi] = unpack_uv(uv);
value = inner_objective(theta, phi, A, lambda_k, ax_tol);
end

function value = inner_objective(theta, phi, A, lambda_k, ax_tol)
x = unit_vector(theta, phi);
Ax = A * x;
gx = norm(Ax, 2);

if gx <= ax_tol
    value = inf;
    return;
end

fx = real(x' * A * x);
value = fx - lambda_k * gx;
end

function x = unit_vector(theta, phi)
x = [cos(theta); exp(1i * phi) * sin(theta)];
end

function [theta, phi] = unpack_uv(uv)
theta = (pi / 4) * (tanh(uv(1)) + 1);
phi = mod(uv(2), 2 * pi);
end

function u = theta_to_u(theta)
scaled = (4 * theta / pi) - 1;
scaled = min(max(scaled, -0.999999999), 0.999999999);
u = atanh(scaled);
end

function value = get_option(opts, field_name, default_value)
if isfield(opts, field_name)
    value = opts.(field_name);
else
    value = default_value;
end
end
