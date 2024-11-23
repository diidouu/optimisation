function [g]=gradient2(R,X,Y,Cx,Cy,sig)
% On stocke les valeurs du vecteur
g=[0,0];
for i=1:length(X)
    Di=sqrt((X(i)-Cx).^2+(Y(i)-Cy).^2);

    g(1)=g(1)-(2*(X(i)-Cx).*(Di-R))/((sig.^2*Di)*(((R-Di)/sig).^2+1));
    g(2)=g(2)-(2*(Y(i)-Cy).*(Di-R))/((sig.^2*Di)*(((R-Di)/sig).^2+1));
end