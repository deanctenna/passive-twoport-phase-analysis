clc; clear; close all;

A = input('Please enter a 2x2 complex matrix A, for example [1 2+1i; -3i 4]: ');
N = 25000;
Z = zeros(N,1);

for k = 1:N
    x = randn(2,1) + 1i*randn(2,1);
    Ax = A * x;
    den = norm(x) * norm(Ax);

    if den > 1e-12
        Z(k) = (x' * A * x) / den;
    else
        Z(k) = NaN;
    end
end

Z = Z(~isnan(Z));

% Only keep the right-half part: Re(Z) >= 0
Z = Z(real(Z) >= 0);

xr = real(Z);
yi = imag(Z);

theta = linspace(0, 2*pi, 500);
xc = cos(theta);
yc = sin(theta);

figure;
scatter(xr, yi, 6, 'filled'); hold on;

% Unit circle boundary: black dashed
plot(xc, yc, 'k--', 'LineWidth', 1.5);

% Outer boundary of sampled points: black solid
if numel(xr) >= 3
    k = boundary(xr, yi, 0.8);
    plot(xr(k), yi(k), 'k-', 'LineWidth', 1.8);
end

xline(0, 'k-', 'LineWidth', 1);
yline(0, 'k-', 'LineWidth', 1);

axis equal;
grid on;
box on;
xlim([-1.1 1.1]);
ylim([-1.1 1.1]);

xlabel('Re(Z)');
ylabel('Im(Z)');
title('Right-half part of Z = (x^* A x) / (||x|| ||Ax||)');
legend('Sample points', 'Unit circle boundary', 'Outer boundary', 'Location', 'best');
