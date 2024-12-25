function [t, y] = run_rk4(f, t0, y0, h, tf)
    % 龙格-库塔方法的实现
    t = t0:h:tf;  % 创建 t 的离散点
    n = length(t);
    y = zeros(length(y0), n);
    y(:, 1) = y0;

    for i = 1:n-1
        k1 = h * f(t(i), y(:, i));
        k2 = h * f(t(i) + h/2, y(:, i) + k1/2);
        k3 = h * f(t(i) + h/2, y(:, i) + k2/2);
        k4 = h * f(t(i) + h, y(:, i) + k3);

        y(:, i+1) = y(:, i) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end
end

