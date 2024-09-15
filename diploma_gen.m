%   Фиксируем размерность задачи

Nq = 10; Np = 10; Nx = 10; Nz = 10;
N = Nx * Nz; M = Nq * Np;

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
y = A * arr_g + w;...моделируем сигнал

%   Данные взяты из Таблицы 1
r = MBIR(y, 2, 2, 1.1, 0.05, fspecial('gaussian', 3, 0.1), 10, 50, M, D, A, 1e-4);
%r = MBIR(y, 2, 2, 1.1, 0.05, fspecial('gaussian', 3, 0.1), 1, 0, M, D, A, 1e-4);
subplot(1, 2, 2);
imshow(to_Matrix(r, Nx, Nz))
colorbar;