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
Cx1(c);Cy1(l);
% Pour l'intervalle [-1,4]x[-1,4]
M2=min(Z2,[],'all');
[l2,c2]=find(Z2 ==M2);
Cx2(c2);Cy2(l2);

figure(3);
hold on;
scatter(donnees.xi,donnees.yi);
contour(Cx2,Cy2,Z2,50);
viscircles([Cx1(c) Cy1(l)],1.5,'Color','r');
viscircles([Cx2(c2) Cy2(l2)],1.5,'Color','b');
title("Solutions approchées pour les intervalles [-1,1]x[-1,2] et [-1,4]x[-1,4]");
axis equal;
hold off;

%% Question 4

g = gradient_CTLS(R,donnees.xi,donnees.yi,0,0);
eps = 10^(-10);
% Comparaison des résultats de la fonction grad et du calcul manuel
g(1)
(fonction_cout(R,donnees.xi,donnees.yi,eps,0) - fonction_cout(R,donnees.xi,donnees.yi,0,0))/eps
g(2)
(fonction_cout(R,donnees.xi,donnees.yi,0,eps) - fonction_cout(R,donnees.xi,donnees.yi,0,0))/eps

%% Question 5

% Création d'une matrice pour contenir les vecteurs du gradient de CTLS en chaque point
n2=30;
Z3=zeros(n2,n2,2);
Cx3=linspace(-1,4,n2);
Cy3=linspace(-1,4,n2);
% Chaque valeur du gradient est calculée et mis dans Z3
for x=1:n2
    for y=1:n2
           g=gradient_CTLS(R,donnees.xi,donnees.yi,Cx3(x),Cy3(y));
           Z3(y,x,1)=g(1);
           Z3(y,x,2)=g(2);
    end
end

figure(4);
hold on;
contour(Cx2,Cy2,Z2,50);   % lignes de niveau de la fonction de coût 
title("Gradient entre [-1,4]x[-1,4]");
grid("on");
colorbar
quiver(Cx3,Cy3,Z3(:,:,1),Z3(:,:,2));  % vecteurs du gradient
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
while(norm(gradient_CTLS(R,donnees.xi,donnees.yi,xk(1),xk(2)))>=esp && i<100)
    evoCTLS=[evoCTLS;fonction_cout(R,donnees.xi,donnees.yi,xk(1),xk(2))];
    evonorm=[evonorm;norm(gradient_CTLS(R,donnees.xi,donnees.yi,xk(1),xk(2)))];
    %Calcul la direction du gradient 
    dir=-gradient_CTLS(R,donnees.xi,donnees.yi,xk(1),xk(2));

    %Calcul alpha grâce à la fonction de Fletcher-Lemarechal
    dist=xk;
    ak=Fletcher(R,xk,dir,donnees.xi,donnees.yi);
    xk=xk+ak*dir;
    i=i+1;
    X=[X;xk];
    evodist=[evodist;norm(dist-xk)];
end

% On affiche les points de la méthode des plus fortes pentes
figure(5);
hold on;
plot(X(:,1),X(:,2),'-');
contour(Cx2,Cy2,Z2,50);
title("Méthode des plus fortes pentes");
hold off;

% On affiche l'évolution de la fonction de coût
figure(6);
subplot(1,3,1);
plot(evoCTLS);
title("Évolution de la fonction de coût")
xlabel("nombre d'itérations"); 
ylabel('valeur fonction cout');

% On affiche l'évolution de la norme du gradient
subplot(1,3,2);
plot(evonorm);
title("Évolution de la norme du gradient");
xlabel("nombre d'itérations");
ylabel('norme du vecteur gradient');

% On affiche l'évolution de la distance entre chaque point 
subplot(1,3,3)
plot(evodist);
title("Évolution de la distance");
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
    donnees = load('measured_points.mat');
    i=0;
    while(norm(gradient_CTLS(R,donnees.xi,donnees.yi,xk(1),xk(2)))>=eps && i<100)
        d=-gradient_CTLS(R,donnees.xi,donnees.yi,xk(1),xk(2));
        ak=Fletcher(R,xk,d,donnees.xi,donnees.yi);
        xk=xk+ak*d;
        i=i+1;
        X=[X;xk];
    end
    plot(X(:,1),X(:,2),'-');
end
contour(Cx2,Cy2,Z2,50);
title("Analyse des itérés avec d'autres points de départ")
hold off;

%% Question 8
%{
% On reprend le xk de la question 6
xk=[0,0]; X=[xk];
i= 0;
I=eye(2); Ak=I;

QuasiNevoCTLS=[]; QuasiNevonorm=[]; QuasiNevodist=[];
while(norm(gradient_CTLS(R,donnees.xi,donnees.yi,xk(1),xk(2)))>=esp && i<100)
    QuasiNevoCTLS=[QuasiNevoCTLS;fonction_cout(R,donnees.xi,donnees.yi,xk(1),xk(2))];
    QuasiNevonorm=[QuasiNevonorm;norm(gradient_CTLS(R,donnees.xi,donnees.yi,xk(1),xk(2)))];
    %Calcul la direction du gradient 
    dir=-Ak*(gradient_CTLS(R,donnees.xi,donnees.yi,xk(1),xk(2))');
  
    %Calcul de x* grâce à l'algorithme de Quasi-Newton
    dist=xk;
    ginit = gradient_CTLS(R,donnees.xi,donnees.yi,xk(1),xk(2));
    ak=Fletcher(R,xk,dir',donnees.xi,donnees.yi);
    xk=xk+ak*dir';
    i=i+1;
    X=[X;xk];
    yk= gradient_CTLS(R,donnees.xi,donnees.yi,xk(1),xk(2)) - ginit;
    dk = (xk-X(k+1,:));
    Hk = (I-(dk*yk')/(dk'*yk))*Ak*(I-(yk*dk')/(dk'*yk))+((dk*dk')/(dk'*yk));
    QuasiNevodist=[QuasiNevodist;norm(dist-xk)];
    end

% On affiche les points de la méthode des plus fortes pentes
figure(13);
hold on;
plot(X(:,1),X(:,2),'-');
contour(Cx2,Cy2,Z2,50);
title("Méthode des plus fortes pentes");
hold off;

% On affiche l'évolution de la fonction de coût
figure(14);
subplot(1,3,1);
plot(QuasiNevoCTLS);
title("Évolution de la fonction de coût")
xlabel("nombre d'itérations"); 
ylabel('valeur fonction cout');

% On affiche l'évolution de la norme du gradient
subplot(1,3,2);
plot(QuasiNevonorm);
title("Évolution de la norme du gradient");
xlabel("nombre d'itérations");
ylabel('norme du vecteur gradient');

% On affiche l'évolution de la distance entre chaque point 
subplot(1,3,3)
plot(QuasiNevodist);
title("Évolution de la distance");
xlabel("nombre d'itérations");
ylabel('distance');
%}

%% Question 9

% On initialise les différents sigma
sig1=10^-3;
sig2=0.1;
sig3=10;

n=100;
Cx5=linspace(-1,4,n);
Cy5=linspace(-1,4,n);

% On  initialise les différentes matrices de la nouvelle fonction de coût
Z5= zeros(n);
Z6= zeros(n);
Z7= zeros(n);

for j=1:n
    for k=1:n
        Z5(j,k)=fonction_cout2(R,donnees.xi,donnees.yi,Cx5(k),Cy5(j),sig1);
        Z6(j,k)=fonction_cout2(R,donnees.xi,donnees.yi,Cx5(k),Cy5(j),sig2);
        Z7(j,k)=fonction_cout2(R,donnees.xi,donnees.yi,Cx5(k),Cy5(j),sig3);
    end
end

% Recherche du minimum pour σ = 0.001
M5=min(Z5,[],'all');
[l5,c5]=find(Z5==M5);
Cx5(c5);
Cy5(l5);

% Contour de la nouvelle fonction de coût pour sigma=10^-3 avec son cercle qui est représenté par le minimum
figure(8);
subplot(1,3,1);
hold on;
contour(Cx5,Cy5,Z5,20);
title("σ = 0.001");
grid("on");
viscircles([Cx5(c5) Cy5(l5)],1.5);
colorbar;
axis equal;

% σ = 0.1
M6=min(Z5,[],'all');
[l6,c6]=find(Z6 ==M6);
Cx5(c6); Cy5(l6);

subplot(1,3,2);
contour(Cx5,Cy5,Z6,20);
title("σ = 0.1");
grid("on");
viscircles([Cx5(c5) Cy5(l5)],1.5);
colorbar;
axis equal;

% σ = 10
M7=min(Z7,[],'all');
[l7,c7]=find(Z7 ==M7);
Cx5(c7);
Cy5(l7);

subplot(1,3,3)
contour(Cx5,Cy5,Z7,20);
title("σ = 10");
grid("on");
viscircles([Cx5(c7) Cy5(l7)],1.5);
hold off;
colorbar;
axis equal;

%% Question 10.(4)

g = gradient_CTLS2(R,donnees.xi,donnees.yi,0,0,sig2);
eps = 10^(-10);

g(1)
(fonction_cout2(R,donnees.xi,donnees.yi,eps,0,sig2) - fonction_cout2(R,donnees.xi,donnees.yi,0,0,sig2))/eps
g(2)
(fonction_cout2(R,donnees.xi,donnees.yi,0,eps,sig2) - fonction_cout2(R,donnees.xi,donnees.yi,0,0,sig2))/eps

%% Question 10.(5)

% Création d'une matrice pour contenir les vecteurs du gradient de CTLS en chaque point
Z8=zeros(n2,n2,2);

for x=1:n2
    for y=1:n2
           % Chaque valeur du gradient est ainsi calculée grâce à gradient_CTLS et est mis dans une des coordonnées du tableau
           g=gradient_CTLS2(R,donnees.xi,donnees.yi,Cx3(x),Cy3(y),sig2);
           Z8(y,x,1)=g(1);
           Z8(y,x,2)=g(2);
    end
end

% Affichage des lignes de niveau de la fonction de coût 
figure(9);
contour(Cx2,Cy2,Z6,50);
title("Gradient de la 2ème fonction CTLS avec σ = 0.1 entre [-1,4]x[-1,4]");
grid("on");
colorbar
hold on;
% Affichage des vecteurs du gradient de la fonction de coût sur le même graphique
quiver(Cx3,Cy3,Z8(:,:,1),Z8(:,:,2));
hold off;
axis equal;

%% Question 10.(6)

xk=[0,0];
X=[xk];
eps=0.000001;
i=0;

% Stockage de chaque évolution dans un tableau pour ensuite les réutiliser pour les afficher
evoCTLS=[];
evonorm=[];
evodist=[];

% Fonction des plus fortes pentes
while(norm(gradient_CTLS2(R,donnees.xi,donnees.yi,xk(1),xk(2),sig2))>=eps && i<100)
    evoCTLS=[evoCTLS;fonction_cout2(R,donnees.xi,donnees.yi,xk(1),xk(2),sig2)];
    evonorm=[evonorm;norm(gradient_CTLS2(R,donnees.xi,donnees.yi,xk(1),xk(2),sig2))];
    dist=xk;
    
    %Calcul la direction du gradient 
    d=-gradient_CTLS2(R,donnees.xi,donnees.yi,xk(1),xk(2),sig2);

    %Calcul alpha grâce à la fonction de Fletcher-Lemarechal
    ak=Fletcher2(R,xk,d,donnees.xi,donnees.yi,sig2);
    xk=xk+ak*d;
    i=i+1;
    X=[X;xk];
    evodist=[evodist;norm(dist-xk)];
end

% On affiche les points de la méthode des plus fortes pentes
figure(10);
hold on;
plot(X(:,1),X(:,2),'-');
contour(Cx5,Cy5,Z6,50);
hold off;
title("Méthode des plus fortes pentes");

% On affiche l'évolution de la fonction de coût
figure(11);
subplot(1,3,1);
plot(evoCTLS);
title("Évolution de la fonction de coût");
xlabel("nombre d'itérations"); 
ylabel('valeur fonction cout');

% On affiche l'évolution de la norme du gradient
subplot(1,3,2);
plot(evonorm);
title("Évolution de la norme du gradient");
xlabel("nombre d'itérations"); 
ylabel('norme du vecteur gradient');

% On affiche l'évolution de la distance entre chaque point 
subplot(1,3,3);
plot(evodist);
title("Évolution de la distance");
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
    eps=0.000001;
    R=1.5;
    donnees = load('measured_points.mat');
    i=0;
    while(norm(gradient_CTLS2(R,donnees.xi,donnees.yi,xk(1),xk(2),sig2))>=eps && i<100)
        d=-gradient_CTLS2(R,donnees.xi,donnees.yi,xk(1),xk(2),sig2);
        ak=Fletcher2(R,xk,d,donnees.xi,donnees.yi,sig2);
        xk=xk+ak*d;
        i=i+1;
        X=[X;xk];
    end
    plot(X(:,1),X(:,2),'-');
end
contour(Cx5,Cy5,Z6,50);
title("Analyse des itérés avec d'autres points de départ");
hold off;
