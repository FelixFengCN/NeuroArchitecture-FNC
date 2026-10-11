function output = triangular_conversion(input, mode, type)
% input : 矩阵或向量
% mode  : 'vec2mat' 或 'mat2vec'
% type  : 'lower' 或 'upper'

% 例：
% output = triangular_conversion(vec, 'vec2mat', 'lower');
% output = triangular_conversion(mat, 'mat2vec', 'upper');

    if strcmp(mode, 'vec2mat')
        % 向量转对称矩阵
        n = solve_n(length(input));
        mat = zeros(n);
        if strcmp(type, 'lower')
            mat(tril(true(n), -1)) = input;
            mat = mat + mat';
        elseif strcmp(type, 'upper')
            mat(triu(true(n), 1)) = input;
            mat = mat + mat';
        else
            error('type must be "lower" or "upper"');
        end
        output = mat;

    elseif strcmp(mode, 'mat2vec')
        % 矩阵转向量
        if ~ismatrix(input) || size(input,1) ~= size(input,2)
            error('Input must be a square matrix');
        end
        n = size(input,1);
        if strcmp(type, 'lower')
            output = input(tril(true(n), -1));
        elseif strcmp(type, 'upper')
            output = input(triu(true(n), 1));
        else
            error('type must be "lower" or "upper"');
        end

    else
        error('mode must be "vec2mat" or "mat2vec"');
    end
end

function n = solve_n(l)
    % 解 n*(n-1)/2 = l，得到矩阵维度 n
    n = (1 + sqrt(1 + 8 * l)) / 2;
    if mod(n,1) ~= 0
        error('Input vector length is invalid: cannot reshape into symmetric matrix');
    end
    n = round(n);
end
