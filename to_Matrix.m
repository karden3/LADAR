function Matrix = to_Matrix (array, numRows, numCols)...превращает матрицу в массив

if(size(array) ~= numRows * numCols)
    error("Error in the size of array's size");
end

Matrix = zeros(numRows, numCols);

for i = 1:numRows
    for j = 1:numCols
        Matrix(i, j) = array((i - 1) * numCols + j);
    end
end
end
