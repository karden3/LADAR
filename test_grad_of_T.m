function test_grad_of_T (Nx, Nz, Nq, Np, SNR, p, alpha1, alpha2, alpha3, eps, bw2)

%   ������ ��������� ���������� ��������� ������
%   �� ������� � ���� "����� �"
%   ������������ ����� �����֣���� ����������
%   ��� ���������� ����� ���������
%
%   g, g1          - ��������� � ��������������� �����������, ��������������
%   Nx, Nz, Nq, Np - ����������� ���������, ������������ ������� ������
%   SNR            - var(y) / var(w), y - ������, w - ���
%   A              - ����������� �������
%   w              - ��������� ���
%   y              - ��������������� ��������� ������
%   p              - ������������ ����������� ����� (\eps ��� N)
%   alpha1         - ��������s ������������� ��� ��

%g = generation_g_T(Nq, Np);...�������������� ����������� � ���� "����� �"
g = bw2;
%g = rotate_Matrix (g, 50, 50, pi/3);    
subplot(2, 5, 1);...�����ģ� � ���� ��� ������ ��������
imshow(g);...��������� ��������� �������
str1 = sprintf('��������� �������');...��������� ���������
title(str1)
colorbar;...��������� �������� �����    

subplot(2, 5, 6);...����� �� ��������� ����
x = linspace(1, 50, 50);...���������� ��� �����
position_of_cut = 3 * floor(Nq / 10) + 1;
plot(x, g(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'green', 'Color', 'green');
xlim([0 50])
ylim([0 2])
title(str1)
grid on

A = generation_A(Nx, Nz, Nq, Np);...�������������� ������� �������

arr_g = to_array(g);...���������� ������� � ������

w = generation_w(A, g, SNR);...���������� ���
y = A * arr_g + w;...���������� ������
%y = zeros(1, 250);
%D_phi = zeros(Nq, Np);
%for ii = 1:100
 %   w = generation_w(A, g, SNR);
  %  
   % phi = generation_phi(Nq, Np);...����� �� ���������� �������������� ���������� ������� 
%
 %   for i = 1:Nq
  %      for k = 1:Np
   %         m = (i - 1) * Np + k;...��� ���� ��������� �� ��������� � ��� ��������� ��� �����
    %        D_phi(m, m) = exp(-1i * phi(i, k));...�� �����������
     %   end
    %end
    
    %y = y + D_phi * A * arr_g + w;
%end

%y = y / 100;    

g2 = (A' * A)^(-1) * A' * y;...��������� ���
g2 = to_Matrix(g2, Nx, Nz);...��������� ������� ������ � ���������� ����
g2 = abs(g2);...��ң� ����������� ��������
var2 = rmsd (arr_g, to_array(g2), Nx * Nz);...������� �������� �������

subplot(2, 5, 2);...�����ģ� � ���� ��� ������ ��������
imshow(g2);...���������� ������� � ���
str2 = sprintf('1/SNR = %1.2f, RMSD = %1.4f', 1/SNR, var2);
title(str2)
colorbar;...��������� �������� �����

subplot(2, 5, 7);
plot(x, g(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'green', 'Color', 'green');
hold on
grid on
xlim([0 50])
ylim([0 2])
plot(x, g2(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'red', 'Color', 'red');
title(str2)
hold off

g3 = (A' * A + alpha1 * eye(size(A' * A)))^(-1) * A' * y;...��� � ��������������
g3 = to_Matrix (g3, Nx, Nz);...���������� ��������� ���
g3 = abs(g3);...��ң� ����������� ��������
var3 = rmsd (arr_g, to_array(g3), Nx * Nz);...������� �������� �������

subplot(2, 5, 3);...�����ģ� � ���� ��� ������� ��������
imshow(g3);...��������� ��������� �������
str3 = sprintf('alpha = %d, RMSD = %1.4f', alpha1, var3);...��������� ���������
title(str3);...����������� ���������
colorbar;...��������� �������� �����   

subplot(2, 5, 8);
plot(x, g(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'green', 'Color', 'green');
hold on
grid on
xlim([0 50])
ylim([0 2])
plot(x, g3(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'red', 'Color', 'red');
title(str3)
hold off

AA = A' * A;...��� ���������� ��������� ������������ ������
g4 = grad_reg_fun_of_lim_var (A, AA, y, p, alpha2, eps, Nx, Nz);...������������� 
g4 = to_Matrix (g4, Nx, Nz);...���������� ��������� ���
g4 = abs(g4);...��ң� ����������� ��������
var4 = rmsd (arr_g, to_array(g4), Nx * Nz);...������� �������� �������

subplot(2, 5, 4);...�����ģ� � ���� ��� ���ף���� ��������
imshow(g4);...��������� ��������� �������
str4 = sprintf('RMSD = %1.4f, alpha = %d, eps = %1.2f', var4, alpha2, eps);...��������� ���������
title(str4);...����������� ���������
colorbar;...��������� �������� �����

subplot(2, 5, 9);
plot(x, g(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'green', 'Color', 'green');
hold on
grid on
xlim([0 50])
ylim([0 2])
plot(x, g4(position_of_cut, :), 'Marker', 'o', 'MarkerFaceColor', 'red', 'Color', 'red');
title(str4)
hold off

g5 = grad_reg_fun_of_lim_var_sign (A, AA, y, p, alpha3, Nx, Nz);...������������� 
g5 = to_Matrix (g5, Nx, Nz);...���������� ��������� ���
g5 = abs(g5);...��ң� ����������� ��������
var5 = rmsd (arr_g, to_array(g5), Nx * Nz);...������� �������� �������

subplot(2, 5, 5);...�����ģ� � ���� ��� ���ף���� ��������
imshow(g5);...��������� ��������� �������
str5 = sprintf('RMSD = %1.4f, alpha = %d', var5, alpha3);...��������� ���������
title(str5);...����������� ���������
colorbar;...��������� �������� �����

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