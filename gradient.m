function grad = gradient(cx, cy, points, R)
    grad = [0, 0];
    for i = 1:length(points.xi)
        Di = sqrt((cx - points.xi(i)).^2 + (cy - points.yi(i)).^2);
        grad(1) = grad(1) + 2*(cx - points.xi(i))*(Di - R)/Di;
        grad(2) = grad(2) + 2*(cy - points.yi(i))*(Di - R)/Di;
    end
end