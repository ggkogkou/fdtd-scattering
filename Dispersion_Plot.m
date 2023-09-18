% Step 1: Define symbolic function f(x)
syms f(x)

% Step 2: Define array 'a'
a = [1/10 1/15 1/20 1/25 1/30];

% Step 3: Loop over 'a' values and plot
hold on; % Enable 'hold' to overlay plots

for i = 1:length(a)
    f(x) = sqrt(2)/pi/a(i)*asin(sqrt(2)^(-1) * sqrt(sin(pi*cos(x)*a(i))^2 ...
        + sin(pi*sin(x)*a(i))^2));
    fplot(f, [0 pi/2]);
end

% Step 4: Customize the plot
title('Plot of numerical(v_p)/real(v_p) for different values');
xlabel('θ_{inc}');
ylabel('numerical(v_p)/real(v_p)');
legend('Δx/λ = 1/10', 'Δx/λ = 1/15', 'Δx/λ = 1/20', 'Δx/λ = 1/25', 'Δx/λ = 1/30');
hold off; % Disable 'hold' after plotting is done
