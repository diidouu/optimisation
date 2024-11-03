function J = cost_function(cx, cy, points, R)
    J = sum((sqrt((points.xi - cx).^2 + (points.yi - cy).^2) - R).^2);
end