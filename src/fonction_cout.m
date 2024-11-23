function[cout]=fonction_cout (R,X, Y,cx,cy)

cout=0;
for i=1:length(X)
    Di=sqrt((X(i)-cx)^2+(Y(i)-cy)^2);
    cout=cout+(Di-R)^2;
end



