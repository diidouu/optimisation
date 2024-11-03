clear all;
close all;
clc;

%% Question 1

points = load('measured_points.mat');

scatter(points.xi, points.yi);

R = 1.5; % Rayon désiré 
n = 20; % nombre d'échantillon

% Etude dans l'intervalle [-1, 1]x[-1, 2]

cx_values = linspace(-1, 1, n);
cy_values = linspace(-1, 2, n);
cz = zeros(n,n);

for i = 1:n
    for j = 1:n
        cx = cx_values(i);
        cy = cy_values(j);
        cz(j,i) = cost_function(cx, cy, points, R);
    end
end

figure;
surf(cx_values, cy_values, cz);
colorbar;
title('Cost function');
xlabel('cx');
ylabel('cy');
axis equal;

[x_min1, y_min1] = min_cout(cz);

% Etude dans l'intervalle [-1, 4]x[-1, 4]

cx_values = linspace(-1, 4, n);
cy_values = linspace(-1, 4, n);
cz = zeros(n,n);

for i = 1:n
    for j = 1:n
        cx = cx_values(i);
        cy = cy_values(j);
        cz(j,i) = cost_function(cx, cy, points, R);
    end
end

figure;
surf(cx_values, cy_values, cz);
colorbar;
title('Cost function');
xlabel('cx');
ylabel('cy');
axis equal;

[x_min2, y_min2] = min_cout(cz);

% les premiers intervalles sont trop petits pour voir les minimas locaux. Cela fonctionne mieux avec les deucièmes intervalles.

%% Question 2

figure;
scatter(points.xi, points.yi);
hold on;
a = viscircles([cx_values(x_min1), cy_values(y_min1)], R, 'Color', 'b');
hold on;
b = viscircles([cx_values(x_min2), cy_values(y_min2)], R, 'Color', 'r');
legend([a, b], '[-1,1]x[-1,2]', '[-1,4]x[-1,4]');
axis([-1, 4, -1, 4]);
axis equal;

%Pour avoir une précision de l'ordre de 10^-4 sur cx, cy et sur R sur l'intervalle donné, 
%on a alors un nombre d'itération de l'ordre de grandeur de (10^4)^3 ce qui est beaucoup trop grand.

%% Question 4

verif_grad(0, 0, points, R)
verif_grad(1, 1, points, R)

%% Question 5

figure;
champs_de_vecteurs = zeros(n,n,2);

for i = 1:n
    for j = 1:n
        cx = cx_values(i);
        cy = cy_values(j);
        grad = gradient(cx, cy, points, R);
        champs_de_vecteurs(j,i,1) = grad(1);
        champs_de_vecteurs(j,i,2) = grad(2);
    end
end

quiver(cx_values, cy_values, champs_de_vecteurs(:,:,1), champs_de_vecteurs(:,:,2));
hold on;
contour(cx_values, cy_values, cz, 50);
axis([-1, 4, -1, 4]);
axis equal;
