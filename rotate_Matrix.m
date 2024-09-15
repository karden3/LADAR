function A = rotate_Matrix (M, numRows, numCols, angle)

%   ������ ��������� ������������ �������
%   ������� �� �������� ����
%   
%   M       - �������� �������
%   numRows - ����� ����� ������� �������
%   numCols - ����� �������� ������� �������
%   angle   - ���� �������� ������ ������� ������� � OXY
%   A       - ���������� �������

R = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];...������ ������� ��������
oldAxes = [0; 0; 1];...��������������� ������ ��� ������ ���������
newAxes = oldAxes;...��������������� ������ ��� ����� ���������
A = zeros(size(M));...�������������� �������� �������
d = [0 0 round(numCols/2); 0 0 round(numRows/2); 0 0 0];...������� ��������

for jj = 1:numRows
    oldAxes(2) = oldAxes(2) + 1;...����� �� �� 
    for ii = 1:numCols
        oldAxes(1) = ii;...����� �� ��
        newAxes = round((eye(3) + d) * R * (eye(3) - d) * oldAxes);...���������� ������� �������� � ���������
        newX = min(max(newAxes(1), 1), 50);...�������������� ������� ��� OX
        newY = min(max(newAxes(2), 1), 50);...�������������� ������� ��� OY
        if (M(oldAxes(1), oldAxes(2)) == 1)
            A(newX, newY) = M(oldAxes(1), oldAxes(2));...���������������
        end
    end
end

%   ��������� ������ ������� ������ �������
for jj = 2:(numRows - 1)
    for ii = 2:(numCols - 1)
        if(A(ii - 1, jj) + A(ii, jj + 1) == 2)
            A(ii, jj) = 1;
        end
    end
end
        

end