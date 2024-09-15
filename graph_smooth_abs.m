function graph_smooth_abs (x_lim, eps, M, N)

%   ������ ��������� ������ ������
%   ����������� ������ � ��������
%   ��� ��������� � ������ \eps
%
%   x_lim - ���������� ����� �� �� �� �������
%   eps   - �������� ������������� ������
%   M, N  - ����������� ������������ ������

x = linspace(-x_lim, x_lim, 1000);...����� �� ��

smooth_ab = zeros(1, 1000);...����� �������� ����������� ������
origin_abs = zeros(1, 1000);...����� �������� �������� ������

for i = 1:1000
    smooth_ab(i) = smooth_abs (x(i), eps, M, N);...������� �������� � i-�� ����
    origin_abs(i) = abs(x(i));
end

plot(x, smooth_ab, 'r');
hold on
grid on
str4 = sprintf('"���������� ������" ��� eps = %1.2f, M = N = %d', eps, M);...��������� ���������
title(str4);...����������� ���������
plot(x, origin_abs, 'g');
hold off

end