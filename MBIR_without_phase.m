Nq = 10; Np = 10; Nx = 10; Nz = 10;
N = Nx * Nz; M = Nq * Np;
gamma = 2;
%   Создаём незашумлённую матрицу прямой задачи
D = generation_D(Nx, Nz, Nq, Np);
g = generation_g_T(Nq, Np);...инициализируем изображение в виде "буквы Т"
    
% Для 10 на 10:
g(5:6, 5:6) = 1;
%g = bw2;
subplot(1, 2, 1);
imshow(g);...изобразим модельное решение
str1 = sprintf('Модельное решение');...формируем заголовок
title(str1)
colorbar;...добавляет цветовую шкалу

A = D;...пробуем без фазы
%A = generation_A(Nx, Nz, Nq, Np);...инициализируем матрицу системы
arr_g = to_array(g);...вытягиваем матрицу в вектор
w = generation_w(A, g, 2);...моделируем шум
w = abs(w);...пробуем вещественнозначный шум
y_0 = A * arr_g;
y = A * arr_g + w;...моделируем сигнал


%y = y_0;
%   Данные взяты из Таблицы 1
r = abs(D' * y / M).^2;...делю на M так как в статье ошибка(?)
sigma_w2 = (M - 1) / M * var(y);...встроенная вариация другая!
sigma_r  = sqrt((M - 1) / M * var(r)) / gamma;

%sigma_w2 = 2;
[r, ~] = iterative_EM_without_phase(y, r, zeros(M, 1), 2, 2, 1.1, 0, 2, A, D, sigma_w2, sigma_r, M, 1e-2);

subplot(1, 2, 2);
imshow(to_Matrix(r, Nx, Nz))
colorbar;