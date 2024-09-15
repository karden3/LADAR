function graph_smooth_abs (x_lim, eps, M, N)

%   Данная программа строит график
%   сглаженного модуля и обычного
%   для сравнения и выбора \eps
%
%   x_lim - предельная точка по ОХ на графике
%   eps   - параметр сглаживающего модуля
%   M, N  - размерность пространства матриц

x = linspace(-x_lim, x_lim, 1000);...сетка по ОХ

smooth_ab = zeros(1, 1000);...сетка значений сглаженного модуля
origin_abs = zeros(1, 1000);...сетка значений обычного модуля

for i = 1:1000
    smooth_ab(i) = smooth_abs (x(i), eps, M, N);...находим значение в i-ом узле
    origin_abs(i) = abs(x(i));
end

plot(x, smooth_ab, 'r');
hold on
grid on
str4 = sprintf('"Сглаженный модуль" при eps = %1.2f, M = N = %d', eps, M);...формируем заголовок
title(str4);...прописываем заголовок
plot(x, origin_abs, 'g');
hold off

end