function RMSD_per_alpha_and_eps (Nx, Nz, Nq, Np, SNR, p, alpha1, alpha2, p_eps1, p_eps2)

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
%   p                  - ограничивает градиентный спуск (N)
%   p_alpha1, p_alpha2 - граничные значения параметра регуляризации
%   p_eps1, p_eps2     - аналогично, граничные значения для \eps

g = generation_g_T(Nq, Np);...инициализируем изображение в виде "буквы Т"
A = generation_A (Nx, Nz, Nq, Np);...инициализируем матрицу системы

arr_g = to_array(g);...вытягиваем матрицу в вектор

w = generation_w (A, g, SNR);...моделируем шум
y = A * arr_g + w;...моделируем сигнал

%p_alpha = p_alpha1 : p_alpha2 : 20;...инициализируем участвующие степени
%p_alpha = 10.^p_alpha;...инициализируем массив параметров регуляризации
p_alpha = linspace(alpha1, alpha2, 20);...линейный масштаб
p_eps   = linspace(p_eps1, p_eps2, 20);
p_eps   = 10.^p_eps;...логарифмический масштаб

RMSD    = zeros(size(p_alpha'*p_eps));...инициализируем массив невязок
AA = A' * A;...для реализации программы градиентного спуска


for i = 1:20
    for j = 1:20
        g1      = grad_reg_fun_of_lim_var (A, AA, y, p, p_alpha(i), p_eps(j), Nx, Nz);...регуляризация
        RMSD(i, j) = rmsd (arr_g, g1, Nx * Nz);
    end
end

%semilogx(p_alpha, RMSD);
surf(p_alpha', p_eps, RMSD);
   
grid on;
str1 = sprintf('1/SNR = %1.2f, %d иттераций', 1/SNR, p);...формируем заголовок
str2 = sprintf('Зависимость нормы разности решений от параметра регуляризации');

ay = gca;...для последующего обращения 
ay.YScale = 'log';...логарифмический масштаб для epsilon

set(gca,'FontSize',15)
xlabel('Параметр регуляризации alpha','FontSize',20)
ylabel('epsilon','FontSize',20)
zlabel('RMSD','FontSize',20)

title   (str1,'FontSize',22)
%subtitle(str2,'FontSize',18)

end