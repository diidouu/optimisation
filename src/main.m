%% Projet d'optimisation continue : Estimation robuste du centre d'une pièce circulaire
%% Auteur : EL GOTE Ismaïl-Ayoub
clear; close all; clc;

%% 
points = load('measured_points.mat');

%% Question 1

n=100; 
R=1.5;
dbtype fonction_cout.m;

% Création de la plage de donnée entre [-1,1]x[-1,2]
cx1=linspace(-1,1,n);
cy1=linspace(-1,2,n);
% Entre [-1,4]x[-1,4]
cx2=linspace(-1,4,n);
cy2=linspace(-1,4,n);

% Calcul de la fonction de coût pour les intervalles [-1, 1] x [-1, 2]
Z1=zeros(n); 
for i = 1:n
    for j = 1:n
        Z1(j, i) = fonction_cout(R,points.xi,points.yi,cx1(i),cy1(j));
    end
end

% Calcul de la fonction de coût pour les intervalles [-1, 4] x [-1, 4]
Z2=zeros(n);
for i = 1:n
    for j = 1:n
        Z2(j, i) = fonction_cout(R,points.xi,points.yi,cx2(i),cy2(j));
    end
end

% Affichage des lignes de niveau du contour de la fonction de coût de [-1,1]x[-1,2]
figure(1);
contour(cx1,cy1,Z1,50);
title("Q1 : Fonction de coût - Intervalle = [-1, 1]x[-1, 2]");
xlabel('cx'); ylabel('cy');
grid("on"); colorbar;

% Intervalle de [-1,4]x[-1,4]
figure(2);
contour(cx2,cy2,Z2,50);
title("Q1 : Fonction de coût - Intervalle = [-1,4]x[-1,4]");
xlabel('cx'); ylabel('cy');
grid("on"); colorbar;

%% Question 2

% Recherche du minimum de la fonction de coût grâce à la matrice
% Pour l'intervalle [-1,1]x[-1,2]

M1=min(Z1,[],'all');
[l,c]=find(Z1 == M1);
cx1(c);cy1(l);
% Pour l'intervalle [-1,4]x[-1,4]
M2=min(Z2,[],'all');
[l2,c2]=find(Z2 ==M2);
cx2(c2);cy2(l2);

figure(3);
hold on;
scatter(points.xi,points.yi);
%contour(cx2,cy2,Z2,50);
viscircles([cx1(c) cy1(l)],1.5,'Color','r');
viscircles([cx2(c2) cy2(l2)],1.5,'Color','g');
title("Solutions approchées pour les intervalles [-1,1]x[-1,2] et [-1,4]x[-1,4]");
legend("rouge = [-1,1]x[-1,2] ; vert = [-1,4]x[-1,4]");
axis equal;
hold off;

%% Question 4

g = gradient(R,points.xi,points.yi,0,0);
eps = 10^(-10);
% Comparaison des résultats de la fonction grad et du calcul manuel
g(1)
(fonction_cout(R,points.xi,points.yi,eps,0) - fonction_cout(R,points.xi,points.yi,0,0))/eps
g(2)
(fonction_cout(R,points.xi,points.yi,0,eps) - fonction_cout(R,points.xi,points.yi,0,0))/eps

%% Question 5

% Création d'une matrice pour contenir les vecteurs du gradient de CTLS en chaque point
n2=30;
Z3=zeros(n2,n2,2);
cx3=linspace(-1,4,n2);
cy3=linspace(-1,4,n2);
% Chaque valeur du gradient est calculée et mis dans Z3
for x=1:n2
    for y=1:n2
           g=gradient(R,points.xi,points.yi,cx3(x),cy3(y));
           Z3(y,x,1)=g(1);
           Z3(y,x,2)=g(2);
    end
end

figure(4);
hold on;
contour(cx2,cy2,Z2,50);   % lignes de niveau de la fonction de coût 
title("Q5 : Gradient entre [-1,4]x[-1,4]");
grid("on");
colorbar
quiver(cx3,cy3,Z3(:,:,1),Z3(:,:,2));  % vecteurs du gradient
hold off;
axis equal;

%% Question 6

%Initialisation de xk pour commencer le programme
xk=[0,0];
X=[xk]; %Stockage de chaque valeur des points dans X

%Recherche à approximer la valeur à 0.000001
esp=0.000001;
i=0;

% Stockage de chaque évolution dans un tableau pour les afficher
evoCTLS=[]; evonorm=[]; evodist=[];

% Fonction des plus fortes pentes
while(norm(gradient(R,points.xi,points.yi,xk(1),xk(2)))>=esp && i<100)
    evoCTLS=[evoCTLS;fonction_cout(R,points.xi,points.yi,xk(1),xk(2))];
    evonorm=[evonorm;norm(gradient(R,points.xi,points.yi,xk(1),xk(2)))];
    %Calcul la direction du gradient 
    dir=-gradient(R,points.xi,points.yi,xk(1),xk(2));

    %Calcul alpha grâce à la fonction de Fletcher-Lemarechal
    dist=xk;
    ak=Fl(R,xk,dir,points.xi,points.yi);
    xk=xk+ak*dir;
    i=i+1;
    X=[X;xk];
    evodist=[evodist;norm(dist-xk)];
end

% On affiche les points de la méthode des plus fortes pentes
figure(5);
hold on;
plot(X(:,1),X(:,2),'-');
contour(cx2,cy2,Z2,50);
title("Q6 : Méthode des plus fortes pentes");
hold off;

% On affiche l'évolution de la fonction de coût
figure(6);
subplot(1,3,1);
plot(evoCTLS);
title("Q6 : Évolution de la fonction de coût")
xlabel("nombre d'itérations"); 
ylabel('valeur fonction cout');

% On affiche l'évolution de la norme du gradient
subplot(1,3,2);
plot(evonorm);
title("Q6 : Évolution de la norme du gradient");
xlabel("nombre d'itérations");
ylabel('norme du vecteur gradient');

% On affiche l'évolution de la distance entre chaque point 
subplot(1,3,3)
plot(evodist);
title("Q6 : Évolution de la distance");
xlabel("nombre d'itérations");
ylabel('distance');

%% Question 7

Xi=[-1,1;-1,0;0,0;1,0;1,2;4,3;4,4];
[l,c]=size(Xi);

figure(7);
hold on;
for k=1:l
    xk=Xi(k,:);
    X=[xk];
    eps=0.000001;
    R=1.5;
    points = load('measured_points.mat');
    i=0;
    while(norm(gradient(R,points.xi,points.yi,xk(1),xk(2)))>=eps && i<100)
        d=-gradient(R,points.xi,points.yi,xk(1),xk(2));
        ak=Fl(R,xk,d,points.xi,points.yi);
        xk=xk+ak*d;
        i=i+1;
        X=[X;xk];
    end
    plot(X(:,1),X(:,2),'-');
end
contour(cx2,cy2,Z2,50);
title("Q7 : Analyse des itérés avec d'autres points de départ")
hold off;

%% Question 8

% On reprend le xk de la question 6
xk=[0,0]; Y=[xk];
i= 0;
I=eye(2); Ak=I;

QuasiNevoCTLS=[]; QuasiNevonorm=[]; QuasiNevodist=[];
while(norm(gradient(R,points.xi,points.yi,xk(1),xk(2)))>=esp && i<100)
    QuasiNevoCTLS=[QuasiNevoCTLS;fonction_cout(R,points.xi,points.yi,xk(1),xk(2))];
    QuasiNevonorm=[QuasiNevonorm;norm(gradient(R,points.xi,points.yi,xk(1),xk(2)))];
    %Calcul la direction du gradient 
    dir=-Ak*(gradient(R,points.xi,points.yi,xk(1),xk(2))');
  
    %Calcul de x* grâce à l'algorithme de Quasi-Newton
    dist=xk;
    ginit = gradient(R,points.xi,points.yi,xk(1),xk(2));
    ak=Fl(R,xk,dir',points.xi,points.yi);
    xk=xk+ak*dir';
    i=i+1;
    Y=[Y;xk];
    yk= gradient(R,points.xi,points.yi,xk(1),xk(2)) - ginit;
    dk = (xk-Y(i,:));
    Hk = (I - (dk' * yk) / (dk' * yk)) * Ak * (I - (yk' * dk) / (dk' * yk)) + ((dk' * dk) / (dk' * yk));
    QuasiNevodist=[QuasiNevodist;norm(dist-xk)];
end

% On affiche les points de la méthode de Qausi-Newton
figure(13);
hold on;
plot(Y(:,1),Y(:,2),'-');
contour(cx2,cy2,Z2,50);
title("Q8 : Méthode de Quasi-Newton");
hold off;

% On affiche l'évolution de la fonction de coût
figure(14);
subplot(1,3,1);
plot(QuasiNevoCTLS);
title("Q8 : Évolution de la fonction de coût")
xlabel("nombre d'itérations"); 
ylabel('valeur fonction cout');

% On affiche l'évolution de la norme du gradient
subplot(1,3,2);
plot(QuasiNevonorm);
title("Q8 : Évolution de la norme du gradient");
xlabel("nombre d'itérations");
ylabel('norme du vecteur gradient');

% On affiche l'évolution de la distance entre chaque point 
subplot(1,3,3)
plot(QuasiNevodist);
title("Q8 : Évolution de la distance");
xlabel("nombre d'itérations");
ylabel('distance');


%% Question 9

% On initialise les différents sigma
sigma1=10^-3;
sigma2=0.1;
sigma3=10;

n=100;
cx5=linspace(-1,4,n);
cy5=linspace(-1,4,n);

% On  initialise les différentes matrices de la nouvelle fonction de coût
Z5= zeros(n);
Z6= zeros(n);
Z7= zeros(n);

for j=1:n
    for k=1:n
        Z5(j,k)=fonction_cout2(R,points.xi,points.yi,cx5(k),cy5(j),sigma1);
        Z6(j,k)=fonction_cout2(R,points.xi,points.yi,cx5(k),cy5(j),sigma2);
        Z7(j,k)=fonction_cout2(R,points.xi,points.yi,cx5(k),cy5(j),sigma3);
    end
end

% Recherche du minimum pour σ = 0.001
M5=min(Z5,[],'all');
[l5,c5]=find(Z5==M5);
cx5(c5);
cy5(l5);

% Contour de la nouvelle fonction de coût pour sigma=10^-3 avec son cercle qui est représenté par le minimum
figure(8);
subplot(1,3,1);
hold on;
contour(cx5,cy5,Z5,20);
title("Q9 : sigma = 0.001");
grid("on");
viscircles([cx5(c5) cy5(l5)],1.5);
colorbar;
axis equal;

% σ = 0.1
M6=min(Z5,[],'all');
[l6,c6]=find(Z6==M6);
cx5(c6); cy5(l6);

subplot(1,3,2);
contour(cx5,cy5,Z6,20);
title("Q9 : sigma = 0.1");
grid("on");
viscircles([cx5(c5) cy5(l5)],1.5);
colorbar;
axis equal;

% σ = 10
M7=min(Z7,[],'all');
[l7,c7]=find(Z7 ==M7);
cx5(c7);
cy5(l7);

subplot(1,3,3)
contour(cx5,cy5,Z7,20);
title("Q9 : sigam = 10");
grid("on");
viscircles([cx5(c7) cy5(l7)],1.5);
hold off;
colorbar;
axis equal;

%% Question 10.(4)

g = gradient2(R,points.xi,points.yi,0,0,sigma2);
epsilon = 10^(-10);

g(1)
(fonction_cout2(R,points.xi,points.yi,eps,0,sigma2) - fonction_cout2(R,points.xi,points.yi,0,0,sigma2))/epsilon
g(2)
(fonction_cout2(R,points.xi,points.yi,0,eps,sigma2) - fonction_cout2(R,points.xi,points.yi,0,0,sigma2))/epsilon

%% Question 10.(5)

% Création d'une matrice pour contenir les vecteurs du gradient de CTLS en chaque point
Z8=zeros(n2,n2,2);

for x=1:n2
    for y=1:n2
           % Chaque valeur du gradient est ainsi calculée grâce à gradient_CTLS et est mis dans une des coordonnées du tableau
           g=gradient2(R,points.xi,points.yi,cx3(x),cy3(y),sigma2);
           Z8(y,x,1)=g(1);
           Z8(y,x,2)=g(2);
    end
end

% Affichage des lignes de niveau de la fonction de coût 
figure(9);
contour(cx2,cy2,Z6,50);
title("Q10.5 : Gradient de la deuxième fonction de coût avec sigma = 0.1 entre [-1,4]x[-1,4]");
grid("on");
colorbar
hold on;
% Affichage des vecteurs du gradient de la fonction de coût sur le même graphique
quiver(cx3,cy3,Z8(:,:,1),Z8(:,:,2));
hold off;
axis equal;

%% Question 10.(6)

xk=[0,0];
X=[xk];
epsilon=0.000001;
i=0;

% Stockage de chaque évolution dans un tableau pour ensuite les réutiliser pour les afficher
evoCTLS=[];
evonorm=[];
evodist=[];

% Fonction des plus fortes pentes
while(norm(gradient2(R,points.xi,points.yi,xk(1),xk(2),sigma2))>=epsilon && i<100)
    evoCTLS=[evoCTLS;fonction_cout2(R,points.xi,points.yi,xk(1),xk(2),sigma2)];
    evonorm=[evonorm;norm(gradient2(R,points.xi,points.yi,xk(1),xk(2),sigma2))];
    dist=xk;
    
    %Calcul la direction du gradient 
    d=-gradient2(R,points.xi,points.yi,xk(1),xk(2),sigma2);

    %Calcul alpha grâce à la fonction de Fletcher-Lemarechal
    ak=Fl2(R,xk,d,points.xi,points.yi,sigma2);
    xk=xk+ak*d;
    i=i+1;
    X=[X;xk];
    evodist=[evodist;norm(dist-xk)];
end

% On affiche les points de la méthode des plus fortes pentes
figure(10);
hold on;
plot(X(:,1),X(:,2),'-');
contour(cx5,cy5,Z6,50);
hold off;
title("Q10.6 : Méthode des plus fortes pentes");

% On affiche l'évolution de la fonction de coût
figure(11);
subplot(1,3,1);
plot(evoCTLS);
title("Q10.6 : Évolution de la fonction de coût");
xlabel("nombre d'itérations"); 
ylabel('valeur fonction cout');

% On affiche l'évolution de la norme du gradient
subplot(1,3,2);
plot(evonorm);
title("Q10.6 : Évolution de la norme du gradient");
xlabel("nombre d'itérations"); 
ylabel('norme du vecteur gradient');

% On affiche l'évolution de la distance entre chaque point 
subplot(1,3,3);
plot(evodist);
title("Q10.6 : Évolution de la distance");
xlabel("nombre d'itérations"); 
ylabel('distance');

%% Question 10.(7)
Xi=[-1,1;-1,0;1,2;1,0]  
[l,c]=size(Xi);
Xi(1,:);
figure(12);
hold on;
for k=1:l
    xk=Xi(k,:);
    X=[xk];
    epsilon=0.000001;
    R=1.5;
    points = load('measured_points.mat');
    i=0;
    while(norm(gradient2(R,points.xi,points.yi,xk(1),xk(2),sigma2))>=epsilon && i<100)
        d=-gradient2(R,points.xi,points.yi,xk(1),xk(2),sigma2);
        ak=Fl2(R,xk,d,points.xi,points.yi,sigma2);
        xk=xk+ak*d;
        i=i+1;
        X=[X;xk];
    end
    plot(X(:,1),X(:,2),'-');
end
contour(cx5,cy5,Z6,50);
title("Q10.7 : Analyse des itérés avec d'autres points de départ");
hold off;
