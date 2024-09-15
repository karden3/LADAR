function array = to_array (Matrix)...построчно матрицу записывает в массив

[numRows, numCols] = size(Matrix);
size_of_array = numRows * numCols;
array = rand(size_of_array, 1);

for i = 1:numRows
    for j = 1:numCols
        array((i - 1) * numCols + j) = Matrix(i, j);
    end
end
end