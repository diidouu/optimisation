function[cout]=fonction_cout2 (R,X,Y,cx,cy,sig)
% Calcul de la 2ème fonction coût à partir de la formule
cout=0;
for i=1:length(X)
    Di=sqrt((X(i)-cx)^2+(Y(i)-cy)^2);
    cout=cout+0.5*log(1+((Di-R)^2)/(sig^2));
end



