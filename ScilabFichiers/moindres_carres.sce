//Acquisition des points
nbPoints = 11;
nbPoints2 = 10;

function res = donnees(n,K)
    for i = linspace(1,n,n)
            res1(i) = 1000+(100*i);
            res2(i) = K(i);
    end
    res = [res1 res2]
endfunction

//a = Acquisitions(i,nbPoints,'oblack');
K = [178,196,214,230,246,261,263,264,265,264,261]
K2 = [193,190,189,188,189,191,193,195,198,201]
a = donnees(nbPoints,K);
disp(a)

//Exercice 3+4
function X = moindresCarres(x,y,n)
    A = [x,ones(size(x,1),1)];
    xT = x;
    for i = (2:1:n)
        xT = xT.*x;
        A = [xT A]
    end
    X = inv(A'*A)*A'*y;
endfunction
function y = afficheCourbe1(a,n)
     plot(a(:,1),a(:,2))
     zoom_rect([400 0 2600 300]);
endfunction
function y = fMC(x,X)
    l = size(X,1)
    y = zeros(1,size(x,2))
    xT = ones(1,size(x,2))
    for i = linspace(l,1,l)
        y = y+xT*X(i);
        xT = xT.*x;
    end
endfunction

//calcul courbe de degré d par les moindres carrés
x = a(:,1);
y = a(:,2);
degre = 2;
X = moindresCarres(x,y,degre);
disp(X)
afficheCourbe1(a,nbPoints)
ech = linspace(0,2500,500);
plot(ech,(X(1)*ech.^2+X(2)*ech+X(3)),'g')

