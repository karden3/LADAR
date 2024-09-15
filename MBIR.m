function r = MBIR(y, gamma, q, p, T, b, N_k, N_l, M, D, A, eps)
%
%   Данная программа иттеративно находит
%   решение, используя EM-алгоритм
%
%   y        - зашумлённые данные
%   gamma    - параметр регуляризации
%   q, p, T  - параметры QGGMRF
%   b        - окрестность пикселей
%   N_k, N_l - число иттераций в циклах
%   M        - размерность задачи
%   D        - незашумлённый оператор
%   A        - оператор с шумом
%
%   Предложенные данные из таблицы 1:
%   gamma = 2
%   q     = 2
%   T     = 0.05
%   b     = fspecial('gaussian', 3, 0.8) (или 0.1?)
%   N_k   = 10
%   N_l   = 300

%   Инициализировать phi через PGA
phi = zeros(M, 1);

for ii = 1:N_l
    r = abs(D_noise(D, phi)' * y / M).^2;
    sigma_w2 = (M - 1) / M * var(y);...встроенная вариация другая!
    sigma_r  = sqrt(sigma_w2) / gamma;
    [~, phi] = iterative_EM(y, r, phi, 2, 2, 1,...
        fspecial('gaussian', 3, 0.8), N_k, A, D, sigma_w2, sigma_r, M, 0);
end
r = abs(D_noise(D, phi)' * y / M).^2;...делю на M так как в статье ошибка(?)
sigma_w2 = (M - 1) / M * var(y);...встроенная вариация другая!
sigma_r  = sqrt((M - 1) / M * var(r)) / gamma;

[r, ~] = iterative_EM(y, r, phi, q, p, T, b, 2, A, D, sigma_w2, sigma_r, M, eps);

    %   Учёт текущей фазы в матрице оператора
    function A = D_noise(D, phi)
        D_phi = zeros(M, M);...первая матрица для А
        phi = to_Matrix(phi, sqrt(M), sqrt(M));
        for i = 1:sqrt(M)
            for k = 1:sqrt(M)
                m = (i - 1) * sqrt(M) + k;...нам надо двигаться по диагонали и вот выражение для этого
                D_phi(m, m) = exp(-1i * abs(phi(i, k)));...по определению
            end
        end
        A = D_phi * D;
    end
end