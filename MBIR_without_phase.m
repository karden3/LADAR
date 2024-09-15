Nq = 10; Np = 10; Nx = 10; Nz = 10;
N = Nx * Nz; M = Nq * Np;
gamma = 2;
%   ������� �������̣���� ������� ������ ������
D = generation_D(Nx, Nz, Nq, Np);
g = generation_g_T(Nq, Np);...�������������� ����������� � ���� "����� �"
    
% ��� 10 �� 10:
g(5:6, 5:6) = 1;
%g = bw2;
subplot(1, 2, 1);
imshow(g);...��������� ��������� �������
str1 = sprintf('��������� �������');...��������� ���������
title(str1)
colorbar;...��������� �������� �����

A = D;...������� ��� ����
%A = generation_A(Nx, Nz, Nq, Np);...�������������� ������� �������
arr_g = to_array(g);...���������� ������� � ������
w = generation_w(A, g, 2);...���������� ���
w = abs(w);...������� ������������������ ���
y_0 = A * arr_g;
y = A * arr_g + w;...���������� ������


%y = y_0;
%   ������ ����� �� ������� 1
r = abs(D' * y / M).^2;...���� �� M ��� ��� � ������ ������(?)
sigma_w2 = (M - 1) / M * var(y);...���������� �������� ������!
sigma_r  = sqrt((M - 1) / M * var(r)) / gamma;

%sigma_w2 = 2;
[r, ~] = iterative_EM_without_phase(y, r, zeros(M, 1), 2, 2, 1.1, 0, 2, A, D, sigma_w2, sigma_r, M, 1e-2);

subplot(1, 2, 2);
imshow(to_Matrix(r, Nx, Nz))
colorbar;