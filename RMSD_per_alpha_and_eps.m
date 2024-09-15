function RMSD_per_alpha_and_eps (Nx, Nz, Nq, Np, SNR, p, alpha1, alpha2, p_eps1, p_eps2)

%   ������ ��������� ��������� ��������� 
%   ����������� RMSD ����������� � ��������� �������
%   �� ��������� \alpha � ��������������� �����
%
%   g, g1              - ��������� � ��������������� �����������, ��������������
%   Nx, Nz, Nq, Np     - ����������� ���������, ������������ ������� ������
%   SNR                - var(y) / var(w), y - ������, w - ���
%   A                  - ����������� �������
%   w                  - ��������� ���
%   y                  - ��������������� ��������� ������
%   p                  - ������������ ����������� ����� (N)
%   p_alpha1, p_alpha2 - ��������� �������� ��������� �������������
%   p_eps1, p_eps2     - ����������, ��������� �������� ��� \eps

g = generation_g_T(Nq, Np);...�������������� ����������� � ���� "����� �"
A = generation_A (Nx, Nz, Nq, Np);...�������������� ������� �������

arr_g = to_array(g);...���������� ������� � ������

w = generation_w (A, g, SNR);...���������� ���
y = A * arr_g + w;...���������� ������

%p_alpha = p_alpha1 : p_alpha2 : 20;...�������������� ����������� �������
%p_alpha = 10.^p_alpha;...�������������� ������ ���������� �������������
p_alpha = linspace(alpha1, alpha2, 20);...�������� �������
p_eps   = linspace(p_eps1, p_eps2, 20);
p_eps   = 10.^p_eps;...��������������� �������

RMSD    = zeros(size(p_alpha'*p_eps));...�������������� ������ �������
AA = A' * A;...��� ���������� ��������� ������������ ������


for i = 1:20
    for j = 1:20
        g1      = grad_reg_fun_of_lim_var (A, AA, y, p, p_alpha(i), p_eps(j), Nx, Nz);...�������������
        RMSD(i, j) = rmsd (arr_g, g1, Nx * Nz);
    end
end

%semilogx(p_alpha, RMSD);
surf(p_alpha', p_eps, RMSD);
   
grid on;
str1 = sprintf('1/SNR = %1.2f, %d ���������', 1/SNR, p);...��������� ���������
str2 = sprintf('����������� ����� �������� ������� �� ��������� �������������');

ay = gca;...��� ������������ ��������� 
ay.YScale = 'log';...��������������� ������� ��� epsilon

set(gca,'FontSize',15)
xlabel('�������� ������������� alpha','FontSize',20)
ylabel('epsilon','FontSize',20)
zlabel('RMSD','FontSize',20)

title   (str1,'FontSize',22)
%subtitle(str2,'FontSize',18)

end