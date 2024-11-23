function [g]=gradient(R,X,Y,Cx,Cy)
% On stocke les valeurs du vecteur
g=[0,0];
for i=1:length(X)
    Di=sqrt((X(i)-Cx)^2+(Y(i)-Cy)^2);
    
    g(1)=g(1)-2*(X(i)-Cx)*(Di-R)/Di;
    g(2)=g(2)-2*(Y(i)-Cy)*(Di-R)/Di;
end