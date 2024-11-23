function [a] = Fl(R,x,d,xi,yi)
% Fonction de Fletcher pour la première fonction CTLS
% On pose les bêtas 1 et 2 ainsi que les alphas i, r et l selon le cours
B1=10^-18;
B2=0.999;
ai=1;
ar=10^18;
al=0;
g=-B1*gradient_CTLS(R,xi,yi,x(1),x(2))*d';

% On établit les conditions de Wolfe
CW1 = fonction_cout(R,xi,yi,x(1)+ai*d(1),x(2)+ai*d(2)) <=fonction_cout(R,xi,yi,x(1),x(2))-ai*g;
CW2 = (gradient_CTLS(R,xi,yi,x(1)+ai*d(1),x(2)+ai*d(2))*d')/(gradient_CTLS(R,xi,yi,x(1),x(2))*d')<= B2;

% On effectue l'algorithme
while (~(CW1 && CW2))
    if(CW1)
        al=ai;
        ai=(al+ar)/2;
    else
        ar=ai;
        ai=(al+ar)/2;
    end
    CW1 = fonction_cout(R,xi,yi,x(1)+ai*d(1),x(2)+ai*d(2)) <= fonction_cout(R,xi,yi,x(1),x(2))-ai*g;
    CW2 = (gradient_CTLS(R,xi,yi,x(1)+ai*d(1),x(2)+ai*d(2))*d')/(gradient_CTLS(R,xi,yi,x(1),x(2))*d') <= B2;
end
a=ai;
end

