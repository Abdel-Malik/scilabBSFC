//Acquisition des points
nbPoints = 11;
nbPtsConso  = 10;
//Points pour moindres carrés
K = [178,196,214,230,246,261,263,264,265,264,261];
ptsConso = [193,190,189,188,189,191,193,195,198,201];

function X = moindresCarres(x,y,n)
    A = [x,ones(size(x,1),1)];
    xT = x;
    for i = (2:1:n)
        xT = xT.*x;
        A = [xT A]
    end
    X = inv(A'*A)*A'*y;
endfunction

function res = donnees(n,K)
    for i = linspace(1,n,n)
            res1(i) = 1000+(100*i);
            res2(i) = K(i);
    end
    res = [res1 res2];
endfunction

function y = afficheCourbe1(a,n)
     plot(a(:,1),a(:,2));
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
matConsoMdreCre = donnees(nbPtsConso,ptsConso);
matPuissMdreCre = donnees(nbPoints,1000*K);

x = matConsoMdreCre(:,1);
x2 = matPuissMdreCre(:,1);
y = matConsoMdreCre(:,2);
y2 = matPuissMdreCre(:,2);
degre = 2;
X = moindresCarres(x,y,degre);
X2 = moindresCarres(x2,y2,degre);
//afficheCourbe1(matPuissMdreCre,nbPoints)
//ech = linspace(0,2500,500);
//couple = ((30/%pi)*fMC(ech,X2))./ech;
//plot(ech,fMC(ech,X2),'g')
//plot(ech,couple,'b')

function res = ff(x,y)
    res = 200+((x/1100-1).^2)+15*(1670/Y)+(x.*(y/1500))/1900;
endfunction

rpm = linspace(0,2500,200)
c = linspace(0,1700,200);
x = rpm;
y = c;

[X,Y] = meshgrid(x,y);
Z = ff(X,Y);
grayplot(x,y,Z);
