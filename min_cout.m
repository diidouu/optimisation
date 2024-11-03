function [x_min, y_min] = min_cout(cz)
    [n, m] = size(cz);
    minimum = cz(1, 1);
    x_min = 1;
    y_min = 1;
    for i = 1:n
        for j = 1:m
            if cz(i, j) < minimum
                minimum = cz(i, j);
                x_min = i;
                y_min = j;
            end
        end
    end
end