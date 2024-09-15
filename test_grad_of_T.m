function test_grad_of_T (Nx, Nz, Nq, Np, SNR, p, alpha1, alpha2, alpha3, eps, bw2)

%   Данная программа моделирует пришедший сигнал
%   от объекта в виде "буквы Т"
%   используется метод сопряжённых градиентов
%   для уменьшения числа иттераций
%
%   g, g1          - модельное и восстановленное изображение, соответственно
%   Nx, Nz, Nq, Np - характерные параметры, определяющие размеры матриц
%   SNR            - var(y) / var(w), y - сигнал, w - шум
%   A              - операторная матрица
%   w              - модельный шум
%   y              - смоделированный пришедший сигнал
%   p              - ограничивает градиентный спуск (\eps или N)
%   alpha1         - параметрs регуляризации для МН

%g = generation_g_T(Nq, Np);...инициализируем изображение в виде "буквы Т"
g = bw2;
%g = rotate_Matrix (g, 50, 50, pi/3);    
subplot(2, 5, 1);...перейдём в поле для первой картинки
imshow(g);...изобразим модельное решение
str1 = sprintf('Модельное решение');...формируем заголовок
title(str1)
colorbar;...добавляет цветовую шкалу    

subplot(2, 5, 6);...здесь мы изобразим срез
x = linspace(1, 50, 50);...необходимо для среза
position_of_cut = 3 * floor(Nq / 10) + 1;
plot(x, g(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'green', 'Color', 'green');
xlim([0 50])
ylim([0 2])
title(str1)
grid on

A = generation_A(Nx, Nz, Nq, Np);...инициализируем матрицу системы

arr_g = to_array(g);...вытягиваем матрицу в вектор

w = generation_w(A, g, SNR);...моделируем шум
y = A * arr_g + w;...моделируем сигнал
%y = zeros(1, 250);
%D_phi = zeros(Nq, Np);
%for ii = 1:100
 %   w = generation_w(A, g, SNR);
  %  
   % phi = generation_phi(Nq, Np);...здесь мы используем предварительно написанную функцию 
%
 %   for i = 1:Nq
  %      for k = 1:Np
   %         m = (i - 1) * Np + k;...нам надо двигаться по диагонали и вот выражение для этого
    %        D_phi(m, m) = exp(-1i * phi(i, k));...по определению
     %   end
    %end
    
    %y = y + D_phi * A * arr_g + w;
%end

%y = y / 100;    

g2 = (A' * A)^(-1) * A' * y;...реализуем МНК
g2 = to_Matrix(g2, Nx, Nz);...переводим искомый вектор к матричному виду
g2 = abs(g2);...берём амплитудное значение
var2 = rmsd (arr_g, to_array(g2), Nx * Nz);...находим вариацию решения

subplot(2, 5, 2);...перейдём в поле для второй картинки
imshow(g2);...изображаем решение с МНК
str2 = sprintf('1/SNR = %1.2f, RMSD = %1.4f', 1/SNR, var2);
title(str2)
colorbar;...добавляем цветовую шкалу

subplot(2, 5, 7);
plot(x, g(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'green', 'Color', 'green');
hold on
grid on
xlim([0 50])
ylim([0 2])
plot(x, g2(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'red', 'Color', 'red');
title(str2)
hold off

g3 = (A' * A + alpha1 * eye(size(A' * A)))^(-1) * A' * y;...МНК с регуляризацией
g3 = to_Matrix (g3, Nx, Nz);...возвращаем матричный вид
g3 = abs(g3);...берём амплитудное значение
var3 = rmsd (arr_g, to_array(g3), Nx * Nz);...находим вариацию решения

subplot(2, 5, 3);...перейдём в поле для третьей картинки
imshow(g3);...изобразим найденное решение
str3 = sprintf('alpha = %d, RMSD = %1.4f', alpha1, var3);...формируем заголовок
title(str3);...прописываем заголовок
colorbar;...добавляет цветовую шкалу   

subplot(2, 5, 8);
plot(x, g(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'green', 'Color', 'green');
hold on
grid on
xlim([0 50])
ylim([0 2])
plot(x, g3(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'red', 'Color', 'red');
title(str3)
hold off

AA = A' * A;...для реализации программы градиентного спуска
g4 = grad_reg_fun_of_lim_var (A, AA, y, p, alpha2, eps, Nx, Nz);...регуляризация 
g4 = to_Matrix (g4, Nx, Nz);...возвращаем матричный вид
g4 = abs(g4);...берём амплитудное значение
var4 = rmsd (arr_g, to_array(g4), Nx * Nz);...находим вариацию решения

subplot(2, 5, 4);...перейдём в поле для четвёртой картинки
imshow(g4);...изобразим найденное решение
str4 = sprintf('RMSD = %1.4f, alpha = %d, eps = %1.2f', var4, alpha2, eps);...формируем заголовок
title(str4);...прописываем заголовок
colorbar;...добавляет цветовую шкалу

subplot(2, 5, 9);
plot(x, g(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'green', 'Color', 'green');
hold on
grid on
xlim([0 50])
ylim([0 2])
plot(x, g4(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'red', 'Color', 'red');
title(str4)
hold off

g5 = grad_reg_fun_of_lim_var_sign (A, AA, y, p, alpha3, Nx, Nz);...регуляризация 
g5 = to_Matrix (g5, Nx, Nz);...возвращаем матричный вид
g5 = abs(g5);...берём амплитудное значение
var5 = rmsd (arr_g, to_array(g5), Nx * Nz);...находим вариацию решения

subplot(2, 5, 5);...перейдём в поле для четвёртой картинки
imshow(g5);...изобразим найденное решение
str5 = sprintf('RMSD = %1.4f, alpha = %d', var5, alpha3);...формируем заголовок
title(str5);...прописываем заголовок
colorbar;...добавляет цветовую шкалу

subplot(2, 5, 10);
plot(x, g(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'green', 'Color', 'green');
hold on
grid on
xlim([0 50])
ylim([0 2])
plot(x, g5(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'red', 'Color', 'red');
title(str5)
hold off

end 