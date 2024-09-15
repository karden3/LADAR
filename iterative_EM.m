function [r, phi] = iterative_EM(y, r0, phi0, q, p, T, b, N_k, A, D, sigma_w2, sigma_r, M, eps)

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

mu0 = r0;...��������� ����������� ��� ���. ��������
r = r0;...��� �������� ����������� ������

psi = phi0; phi = psi;

flag = 0;...����� ������ �� eps ��� N_k
if eps > 0
    flag = 1;
end
    
for ii = 1:N_k
    %   ��������� ������� ����������� ������� ����� ������������ (38)
    %!  �������� ��� ������������� ��������� - z
    %DD = D_noise(D, phi);
    
    
    DD = D;...������� ��� ����
    
    %D_diag_r = D_diag(r);
    %p_theta = @(g)(norm(y - DD * g) ^ 2 / sigma_w2 + g' * D_diag_r ^ (-1) * g) / 1e16;
    %[~, mu] = fminsearch(p_theta, mu0);
    %mu0 = mu;
    
    C = C_fun(sigma_w2, r0);...������� ����������
    mu = C * DD * y / sigma_w2;
    mu0 = mu;
    
    for jj = 1:M
%         if jj == 1
%             r_0 = 0;
%         else
%             r_0 = r(jj - 1);
%         end
        f_0 = @(r_j)f(r_j, r, C, mu0, jj, b);
        
        r(jj) = fminbnd(f_0, 0, 1);
        %[~, r(jj)] = fminsearch(f_0, r(jj));
        
        r(jj) = abs(r(jj));
    end
    
    sigma_w2 = 1 / M * y' * y - 2 / M * real(y' * DD * mu0);
    for jj = 1:M
        sigma_w2 = sigma_w2 + C(jj, jj) + abs(mu0(jj)) ^ 2;
    end
    sigma_w2 = abs(sigma_w2);
    
    %   ������� p = \overline{1, N} - ��� ������� �����������
    %   �������� ���� ����
%     for jj = 1:M
%         f_psi_sample = @(psi_j) f_psi(psi_j, psi, y, jj, D, mu);
%         psi(jj) = fminbnd(f_psi_sample, 0, 2 * pi);
%         %phi(jj) = exp(psi(jj));     
%         phi(jj) = exp(-1j * psi(jj));
%     end
    
    %   � ������ ��������� �� eps
    if flag == 1
        eps_cur = norm(r - r0) / norm(r0);
        if eps_cur > eps
            N_k = N_k + 1;
        else
            flag = 0;
        end
    end
    r0 = r;
end

%   ���� ���������������� ������� ��� argmin





    %   ��������� �������, ������� ������� ������������� r_{s} (43)
    function f_r_s = f(r_j, r, C, mu, s, b)
        %f_r_s = log(abs(r_j)) + (C(s, s) + abs(mu(s)) ^ 2) / abs(r_j);
        
        f_r_s = log(r_j) + (C(s, s) + abs(mu(s)) ^ 2) / r_j;
        %   ������� \rho ������ � (36)
        rho = @(delta) abs(delta) ^ p / p * ((delta / T) ^ (q - p)) / (1 + (delta / T) ^ (q - p));
        [dim, ~] = size(r);
        dim = sqrt(dim);
        r_new = zeros(dim + 2);...��� ������������� ������� ����������
        r_new(2:dim + 1, 2:dim + 1) = to_Matrix(r, dim, dim);...��������� ������ - ������
        % (i - 1) * numCols + j
        pos_j = mod(s, dim); 
        if pos_j == 0
            pos_j = dim;
        end
        pos_i = (s - pos_j) / dim + 1;
%         r_new = to_Matrix(r, dim, dim);
%         for i1 = 1:dim
%             for j1 = 1:dim
%                 f_r_s = f_r_s + b(mod(pos_i - i1), mod(pos_j - j_1)) *...
%                     rho(abs(r_new(i1, j1) - r_j)) / sigma_r ^ p;
%             end
%         end
        for i1 = 1:3
            for j1 = 1:3
                f_r_s = f_r_s + b(i1, j1) *...
                    rho(abs(r_new(pos_i - 1 + i1, pos_j - 1 + j1) - r_j)) / sigma_r ^ p;
            end
        end
%         for ss = 1:M
%             f_r_s = f_r_s + b(s, ss) * rho(abs((r(ss) - r_j)) / sigma_r ^ p);
%         end
    end

    %   ��������� ������� ��� ���ޣ�� argmin ���� (47)
    function res_psi = f_psi(psi_j, psi, y, p_sample, D, mu)
        mat = D * mu;
        % ��� �� ����������� ����������????
        res_psi = - real(y(p_sample, :)' * exp(- 1j * psi_j) * mat(p_sample, :));
        res_psi = abs(res_psi);
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