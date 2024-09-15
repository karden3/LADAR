function RMSD_per_itteration (Nx, Nz, Nq, Np, SNR1, SNR2, N, alpha, eps)

%   ������ ��������� ������ ������ ������� �����������
%   � ���������� ������� �� ���������� ���������
%   � ����������� ������, ��������������� �������������
%   �� ��������� ������� ������������ ��������
%
%   g      - ����������������� �����������
%   �A     - ������� ���������, ���������� �� ����������������� A' * A
%   y      - ��������� ������
%   g0     - ��������� ������
%   N      - ���������� ���������
%   alpha  - �������� ������������� ��������
%   ro     - ������� �� ������� ���������
%   eps    - �������� ������������� ������
%   Nx, Nz - ������� ������� g
%   w      - ��������������� ������ - �������� �������������
%   RMSD   - ������ ������� �������� ������� �� ������ ���������

RMSD1 = zeros(1, N); 
RMSD2 = zeros(1, N); 

%      ���������� ������
g = generation_g_T(Nq, Np);...�������������� ����������� � ���� "����� �"
A = generation_A (Nx, Nz, Nq, Np);...�������������� ������� �������
arr_g = to_array(g);...���������� ������� � ������
w = generation_w (A, g, SNR1);...���������� ���
y = A * arr_g + w;...���������� ������

%      ��������� ����������� �����
AA = A' * A * 2;...��� �������� � ��ԣ��� (�� ����������� 2)
p = zeros(Nx * Nz, 1);
g = zeros(Nx * Nz, 1);...��������� ����������� - �������
w = delta_omega(g, eps, Nx, Nz);

for s = 1:1:N
    if s == 1
        r = AA * g - 2 * A' * y + alpha * w;...2 �� �����������
    else
        r = r - q./ dot(p,q);...dot(a, b) == a' * b
    end
    p = p + r./ dot(r,r);
    w = delta_omega (p, eps, Nx, Nz);
    q = AA * p + alpha * w;
    g = g - p./ dot(p,q);
    RMSD1(s) = rmsd (arr_g, to_array(g), Nx * Nz);...������� �������� �������
end

%   ���������� � ������ SNR
w = generation_w (A, g, SNR2);...���������� ���
y = A * arr_g + w;...���������� ������

%      ��������� ����������� �����
p = zeros(Nx * Nz, 1);
g = zeros(Nx * Nz, 1);...��������� ����������� - �������
w = delta_omega (g, eps, Nx, Nz);

for s = 1:1:N
    if s == 1
        r = AA * g - 2 * A' * y + alpha * w;...2 �� �����������
    else
        r = r - q./ dot(p,q);...dot(a, b) == a' * b
    end
    p = p + r./ dot(r,r);
    w = delta_omega (p, eps, Nx, Nz);
    q = AA * p + alpha * w;
    g = g - p./ dot(p,q);
    RMSD2(s) = rmsd (arr_g, to_array(g), Nx * Nz);...������� �������� �������
end

%      ������ ������ ��������
s = 1:N;
plot(s, RMSD1, 'b');
hold on
grid on
plot(s, RMSD2, 'g');
str4 = sprintf('RMSD(s), eps = %1.5f, alpha = %d', eps, alpha);...��������� ���������
title(str4);...����������� ���������
str1 = sprintf('1/SNR = %1.4f', 1/SNR1);
str2 = sprintf('1/SNR = %1.4f', 1/SNR2);
legend(str1, str2)
ylim([0 0.3])
hold off

end