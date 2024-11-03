function diff = verif_grad(x0, y0, points, R)
    grad = gradient(x0, y0, points, R);
    h = 1e-6;
    a = (cost_function(x0 + h, y0, points, R) - cost_function(x0, y0, points, R))/h;
    b = (cost_function(x0, y0 + h, points, R) - cost_function(x0, y0, points, R))/h;
    diff = abs([a, b] - grad);
end