function [r, phi] = iterative_EM_without_phase(y, r0, phi0, q, p, T, b, N_k, A, D, sigma_w2, sigma_r, M, eps)

%   ������ ������� ������������ ���� 
%   EM-��������� ��� MAP ������ r (��. fig. 6)
%
%   M        - ����������� ������
%   y        - ���������� �����̣���� ������ (y \in C ^ M)
%   r0       - '����֣������', ���������� � �������� ���������� 
%              ����������� ��� ���������� ��������� (r0 \in R ^ N)
%   phi0     - ������� ������ ��-�� ������� ���������, ����������
%              � �������� ���������� ����������� ��� ��
%              ���������� ��������� (\phi0 \in R ^ N)
%   sigma_r  - ��������� ����֣������, ������� � ������
%              ������ ��������� (��. (50)) (sigma_r \in R)
%   sigma_w2 - ������� ��������� �����̣���� ������ - y,
%              ���������� ����������� ����
%   q, T, b  - QGGMRF - ��������� (��. (36) � ��� (51))
%   N_k      - ����� ��������� �� ���������� �������� ��������
%   r, phi   - ���������� ������ ��� ����֣������ � ������� ������
%   gamma    - ����������� ��������� ������������� ��������� (��. (50))
%   A        - ������� ��������� ��������� ������ ������
%   D        - �������̣���� ������� ��������� ��������� ������ ������ A
%   eps      - �������� ���������
%
%   �������������� ����� ��������� � r �������� (50)
%   ������ ������, �� ������� �����

%   r0 = abs(A' * y);...!��������
%   gamma = ...;
%   sigma_w2 = (M - 1) / M * var(y);
%   sigma_r  = sqrt(sigma_w2) / gamma;

r = r0;...��� �������� ����������� ������

psi = phi0; phi = psi;
eps_cur = 2 * eps;

while eps_cur > eps
    
    DD = D;...������� ��� ����;
    
    C = C_fun(sigma_w2, r0);...������� ����������
    mu = C * DD * y / sigma_w2;
    mu0 = mu;
    
    for jj = 1:M
        f_0 = @(r_j)f(r_j, r, C, mu0, jj, b);
        
        r(jj) = fminbnd(f_0, 0, 1);
        %[~, r(jj)] = fminsearch(f_0, r0(jj));
        
        r(jj) = abs(r(jj));
    end
    
    sigma_w2 = 1 / M * y' * y - 2 / M * real(y' * DD * mu0);
    for jj = 1:M
        sigma_w2 = sigma_w2 + C(jj, jj) + abs(mu0(jj)) ^ 2;
    end
    %sigma_w2 = abs(sigma_w2);
    eps_cur = norm(r - r0) / norm(r0);
    r0 = r;
end

%   ���� ���������������� ������� ��� argmin





    %   ��������� �������, ������� ������� ������������� r_{s} (43)
    function f_r_s = f(r_j, r, C, mu, s, b)
        %f_r_s = log(abs(r_j)) + (C(s, s) + abs(mu(s)) ^ 2) / abs(r_j);
        
        f_r_s = log(r_j) + (C(s, s) + abs(mu(s)) ^ 2) / r_j;
        %   ������� \rho ������ � (36)
        %rho = @(delta) abs(delta) ^ p / p * ((delta / T) ^ (q - p)) / (1 + (delta / T) ^ (q - p));
        %for ss = 1:M
        %    f_r_s = f_r_s + b(s, ss) * rho(abs((r(ss) - r_j)) / sigma_r);
        %end
    end

    %   ��������������� ������� g
    function D_res = D_diag(r)
        [dim, ~] = size(r);
        D_res = zeros(dim);
        for kk = 1:dim
            D_res(kk, kk) = r(kk); 
        end
    end
    
    %   ������ ������� ���������� C (41)
    function C_res = C_fun(sigma_w2, r)
        [dim, ~] = size(r);
        C_res = zeros(dim);
        for kk = 1:dim
            C_res(kk, kk) = 1 / (dim + sigma_w2 / r(kk));
        end
        C_res = C_res * sigma_w2;
    end


end