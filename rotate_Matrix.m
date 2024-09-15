function A = rotate_Matrix (M, numRows, numCols, angle)

%   Данная программа осуществляет поворот
%   матрицы на заданный угол
%   
%   M       - исходная матрица
%   numRows - число строк входной матрицы
%   numCols - число столбцрв входной матрицы
%   angle   - угол поворота против часовой стрелки в OXY
%   A       - полученная матрица

R = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];...задаём матрицу поворота
oldAxes = [0; 0; 1];...вспомогательный вектор для старых координат
newAxes = oldAxes;...вспомогательный вектор для новых координат
A = zeros(size(M));...инициализируем выходную матрицу
d = [0 0 round(numCols/2); 0 0 round(numRows/2); 0 0 0];...матрица переноса

for jj = 1:numRows
    oldAxes(2) = oldAxes(2) + 1;...сдвиг по ОУ 
    for ii = 1:numCols
        oldAxes(1) = ii;...сдвиг по ОХ
        newAxes = round((eye(3) + d) * R * (eye(3) - d) * oldAxes);...используем матрицу поворота и округляем
        newX = min(max(newAxes(1), 1), 50);...переопределяем текущую ось OX
        newY = min(max(newAxes(2), 1), 50);...переопределяем текущую ось OY
        if (M(oldAxes(1), oldAxes(2)) == 1)
            A(newX, newY) = M(oldAxes(1), oldAxes(2));...переопределение
        end
    end
end

%   заполняем пустые пиксели внутри объекта
for jj = 2:(numRows - 1)
    for ii = 2:(numCols - 1)
        if(A(ii - 1, jj) + A(ii, jj + 1) == 2)
            A(ii, jj) = 1;
        end
    end
end
        

end