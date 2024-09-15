%   ��������� ����������� ������

Nq = 10; Np = 10; Nx = 10; Nz = 10;
N = Nx * Nz; M = Nq * Np;

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
y = A * arr_g + w;...���������� ������

%   ������ ����� �� ������� 1
r = MBIR(y, 2, 2, 1.1, 0.05, fspecial('gaussian', 3, 0.1), 10, 50, M, D, A, 1e-4);
%r = MBIR(y, 2, 2, 1.1, 0.05, fspecial('gaussian', 3, 0.1), 1, 0, M, D, A, 1e-4);
subplot(1, 2, 2);
imshow(to_Matrix(r, Nx, Nz))
colorbar;