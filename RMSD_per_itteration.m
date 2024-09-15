function RMSD_per_itteration (Nx, Nz, Nq, Np, SNR1, SNR2, N, alpha, eps)

%   Данная программа строит график невязки полученного
%   и модельного решений от количества иттераций
%   в градиентном спуске, осуществляющего регуляризацию
%   на множестве функций ограниченной вариации
%
%   g      - восстанавливаемое изображение
%   АA     - матрица оператора, умноженная на транспонированную A' * A
%   y      - модельный сигнал
%   g0     - начальные данные
%   N      - количество иттераций
%   alpha  - параметр регуляризации Тихонова
%   ro     - невязка на текущей иттерации
%   eps    - параметр сглаживающего модуля
%   Nx, Nz - размеры матрицы g
%   w      - вспомогательный вектор - градиент регуляризации
%   RMSD   - массив искомых значений невязки от номера иттерации

RMSD1 = zeros(1, N); 
RMSD2 = zeros(1, N); 

%      Моделируем данные
g = generation_g_T(Nq, Np);...инициализируем изображение в виде "буквы Т"
A = generation_A (Nx, Nz, Nq, Np);...инициализируем матрицу системы
arr_g = to_array(g);...вытягиваем матрицу в вектор
w = generation_w (A, g, SNR1);...моделируем шум
y = A * arr_g + w;...моделируем сигнал

%      Реализуем градиентный спуск
AA = A' * A * 2;...для сходства с Артёмом (из производной 2)
p = zeros(Nx * Nz, 1);
g = zeros(Nx * Nz, 1);...начальное приближение - нулевое
w = delta_omega(g, eps, Nx, Nz);

for s = 1:1:N
    if s == 1
        r = AA * g - 2 * A' * y + alpha * w;...2 из производной
    else
        r = r - q./ dot(p,q);...dot(a, b) == a' * b
    end
    p = p + r./ dot(r,r);
    w = delta_omega (p, eps, Nx, Nz);
    q = AA * p + alpha * w;
    g = g - p./ dot(p,q);
    RMSD1(s) = rmsd (arr_g, to_array(g), Nx * Nz);...находим вариацию решения
end

%   Аналогично с другим SNR
w = generation_w (A, g, SNR2);...моделируем шум
y = A * arr_g + w;...моделируем сигнал

%      Реализуем градиентный спуск
p = zeros(Nx * Nz, 1);
g = zeros(Nx * Nz, 1);...начальное приближение - нулевое
w = delta_omega (g, eps, Nx, Nz);

for s = 1:1:N
    if s == 1
        r = AA * g - 2 * A' * y + alpha * w;...2 из производной
    else
        r = r - q./ dot(p,q);...dot(a, b) == a' * b
    end
    p = p + r./ dot(r,r);
    w = delta_omega (p, eps, Nx, Nz);
    q = AA * p + alpha * w;
    g = g - p./ dot(p,q);
    RMSD2(s) = rmsd (arr_g, to_array(g), Nx * Nz);...находим вариацию решения
end

%      Строим график вариации
s = 1:N;
plot(s, RMSD1, 'b');
hold on
grid on
plot(s, RMSD2, 'g');
str4 = sprintf('RMSD(s), eps = %1.5f, alpha = %d', eps, alpha);...формируем заголовок
title(str4);...прописываем заголовок
str1 = sprintf('1/SNR = %1.4f', 1/SNR1);
str2 = sprintf('1/SNR = %1.4f', 1/SNR2);
legend(str1, str2)
ylim([0 0.3])
hold off

end