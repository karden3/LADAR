function RMSD_per_alpha (Nx, Nz, Nq, Np, SNR, p, p_alpha1, p_alpha2, eps)

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
%   p                  - ������������ ����������� ����� (\eps ��� N)
%   p_alpha1, p_alpha2 - ��������� �������� ������� ���������
%                        ������������� � ��������������� �����

g = generation_g_T(Nq, Np);...�������������� ����������� � ���� "����� �"
A = generation_A (Nx, Nz, Nq, Np);...�������������� ������� �������

arr_g = to_array(g);...���������� ������� � ������

w = generation_w (A, g, SNR);...���������� ���
y = A * arr_g + w;...���������� ������

%p_alpha = p_alpha1 : p_alpha2 : 20;...�������������� ����������� �������
%p_alpha = 10.^p_alpha;...�������������� ������ ���������� �������������
p_alpha = linspace(p_alpha1, p_alpha2, 20);

RMSD    = zeros(size(p_alpha));...�������������� ������ �������
AA = A' * A;...��� ���������� ��������� ������������ ������

for i = 1:20
    g1      = grad_reg_fun_of_lim_var (A, AA, y, p, p_alpha(i), eps, Nx, Nz);...�������������
    RMSD(i) = rmsd (arr_g, g1, Nx * Nz);
end

%semilogx(p_alpha, RMSD);
plot(p_alpha, RMSD);

grid on;
str1 = sprintf('epsilon = %1.4f, 1/SNR = %1.2f, %d ���������', eps, 1/SNR, p);...��������� ���������
str2 = sprintf('����������� ����� �������� ������� �� ��������� �������������');
set(gca,'FontSize',15)
xlabel('�������� ������������� alpha','FontSize',20)
ylabel('RMSD','FontSize',20)

title   (str1,'FontSize',22)
subtitle(str2,'FontSize',18)


end