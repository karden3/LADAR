function RMSD_per_alpha (Nx, Nz, Nq, Np, SNR, p, p_alpha1, p_alpha2, eps)

%   Данная программа позволяет построить 
%   зависимость RMSD полученного и модельных решений
%   от параметра \alpha в логарифмической шкале
%
%   g, g1              - модельное и восстановленное изображение, соответственно
%   Nx, Nz, Nq, Np     - характерные параметры, определяющие размеры матриц
%   SNR                - var(y) / var(w), y - сигнал, w - шум
%   A                  - операторная матрица
%   w                  - модельный шум
%   y                  - смоделированный пришедший сигнал
%   p                  - ограничивает градиентный спуск (\eps или N)
%   p_alpha1, p_alpha2 - граничные значения степени параметра
%                        регуляризации в логарифмической шкале

g = generation_g_T(Nq, Np);...инициализируем изображение в виде "буквы Т"
A = generation_A (Nx, Nz, Nq, Np);...инициализируем матрицу системы

arr_g = to_array(g);...вытягиваем матрицу в вектор

w = generation_w (A, g, SNR);...моделируем шум
y = A * arr_g + w;...моделируем сигнал

%p_alpha = p_alpha1 : p_alpha2 : 20;...инициализируем участвующие степени
%p_alpha = 10.^p_alpha;...инициализируем массив параметров регуляризации
p_alpha = linspace(p_alpha1, p_alpha2, 20);

RMSD    = zeros(size(p_alpha));...инициализируем массив невязок
AA = A' * A;...для реализации программы градиентного спуска

for i = 1:20
    g1      = grad_reg_fun_of_lim_var (A, AA, y, p, p_alpha(i), eps, Nx, Nz);...регуляризация
    RMSD(i) = rmsd (arr_g, g1, Nx * Nz);
end

%semilogx(p_alpha, RMSD);
plot(p_alpha, RMSD);

grid on;
str1 = sprintf('epsilon = %1.4f, 1/SNR = %1.2f, %d иттераций', eps, 1/SNR, p);...формируем заголовок
str2 = sprintf('Зависимость нормы разности решений от параметра регуляризации');
set(gca,'FontSize',15)
xlabel('Параметр регуляризации alpha','FontSize',20)
ylabel('RMSD','FontSize',20)

title   (str1,'FontSize',22)
subtitle(str2,'FontSize',18)


end