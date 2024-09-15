function r = MBIR(y, gamma, q, p, T, b, N_k, N_l, M, D, A, eps)
%
%   ������ ��������� ����������� �������
%   �������, ��������� EM-��������
%
%   y        - �����̣���� ������
%   gamma    - �������� �������������
%   q, p, T  - ��������� QGGMRF
%   b        - ����������� ��������
%   N_k, N_l - ����� ��������� � ������
%   M        - ����������� ������
%   D        - �������̣���� ��������
%   A        - �������� � �����
%
%   ������������ ������ �� ������� 1:
%   gamma = 2
%   q     = 2
%   T     = 0.05
%   b     = fspecial('gaussian', 3, 0.8) (��� 0.1?)
%   N_k   = 10
%   N_l   = 300

%   ���������������� phi ����� PGA
phi = zeros(M, 1);

for ii = 1:N_l
    r = abs(D_noise(D, phi)' * y / M).^2;
    sigma_w2 = (M - 1) / M * var(y);...���������� �������� ������!
    sigma_r  = sqrt(sigma_w2) / gamma;
    [~, phi] = iterative_EM(y, r, phi, 2, 2, 1,...
        fspecial('gaussian', 3, 0.8), N_k, A, D, sigma_w2, sigma_r, M, 0);
end
r = abs(D_noise(D, phi)' * y / M).^2;...���� �� M ��� ��� � ������ ������(?)
sigma_w2 = (M - 1) / M * var(y);...���������� �������� ������!
sigma_r  = sqrt((M - 1) / M * var(r)) / gamma;

[r, ~] = iterative_EM(y, r, phi, q, p, T, b, 2, A, D, sigma_w2, sigma_r, M, eps);

    %   �ޣ� ������� ���� � ������� ���������
    function A = D_noise(D, phi)
        D_phi = zeros(M, M);...������ ������� ��� �
        phi = to_Matrix(phi, sqrt(M), sqrt(M));
        for i = 1:sqrt(M)
            for k = 1:sqrt(M)
                m = (i - 1) * sqrt(M) + k;...��� ���� ��������� �� ��������� � ��� ��������� ��� �����
                D_phi(m, m) = exp(-1i * abs(phi(i, k)));...�� �����������
            end
        end
        A = D_phi * D;
    end
end