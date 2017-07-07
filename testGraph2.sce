//Acquisition des points
nbPoints = 11;
nbPtsConso  = 10;
//Points pour moindres carrés
ptsPuiss = [178,196,214,230,246,261,263,264,265,264,261];
ptsConso = [193,190,189,188,189,191,193,195,198,201];

//Calcul la solution au sens des moindres carrés
//x : un vecteur de points (données)
//y : un vecteur de points (valeur) ax^n+bx^(n-1)+..p = y
//n : l'ordre du modèle
function X = moindresCarres(x,val,ordre)
    A = [x,ones(size(x,1),1)];
    xT = x;
    for i = (2:1:ordre)
        xT = xT.*x;
        A = [xT A]
    end
    X = inv(A'*A)*A'*val;
endfunction

//création une matrice [x y] de n éléments ; pour le calcul des moindres carrés
//démarre à 1000, echantillonage tout les 100. 
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
matPuissMdreCre = donnees(nbPoints,ptsPuiss);

x = matConsoMdreCre(:,1);
x2 = matPuissMdreCre(:,1);
y = matConsoMdreCre(:,2);
y2 = matPuissMdreCre(:,2);
degre = 2;
X = moindresCarres(x,y,degre);
X2 = moindresCarres(x2,y2,degre);
p = 100
//afficheCourbe1(matPuissMdreCre,nbPoints)
ech = linspace(0,2500,p);
couple = ((30/%pi)*fMC(ech,X2))./ech;
puissance = fMC(ech,X2);


//--------------------------------
function res = n(x)
    res = x ./max(x);
endfunction
function res = ff(x,y,A)
    res = (A(1).*(x.^2))+(A(2).*(y))+A(3);
endfunction

function res = matritialise(a,i)
    res = []
    for i = (1:1:i)
        res = [res a']
    end
endfunction
function res = gg(x,y,A,alpha,a,mcP)
    res = ((alpha+1)*ff(x,y,A)) -matritialise(((2*alpha.*y)./((mcP(1)*(a.^2))+(mcP(2)*a)+mcP(3))),size(ff(x,y,A),1)).*ff(x,y,A) + matritialise(((2*alpha.*y.^2)./(((mcP(1)*(a.^2))+(mcP(2)*a)+mcP(3)).^2)),size(ff(x,y,A),1)).*ff(x,y,A);
endfunction

function res = matriceVal3D(x,mcP,mcC)
    res = (x.^2)';
    p = ((mcP(1)*(x.^2))+(mcP(2)*x)+mcP(3));
    p = p';
    c = ((mcC(1)*(x.^2))+(mcC(2)*x)+mcC(3));
    c=c';
    res = [res (p) ones(size(x,2),1)];
    res = inv(res'*res)*res'*c;
endfunction

function res = afficheConsoPC(x,y,Z,X)
    res = 0;
    for i = (1:1:p)
        a = fMC(x(i),X);
            q = 0;
        for j = (1:2:p)
            if((Z(i,j) >= a-1.5) & (Z(i,j) <= a+1.5) & q < 2) then
                plot(x(i),y(j),'x');
                res = res + 1;
                q = q+1;
            end
        end
    end
endfunction
rpm = linspace(0,2500,p)
c = linspace(1,1700,p);
x = rpm;
y = X2(1)*x.^2+X2(2)*x+X2(3);
a = linspace(600,2500,6);

[A,B] = meshgrid(x,y);
Z = ff(A,B,matriceVal3D(a,X2,X));
f=gcf();f.color_map=hotcolormap(32);xtitle("Graphique d interpolation d un BSFC diesel : f(x,y)=ax²+by+c","regime moteur (tr/min)","puissance fourni (kW)")
zoom_rect([0 0 2600 300]);
colorbar(min(Z'),(max(Z')));
grayplot(x,y,Z');
plot(rpm,puissance)
r = afficheConsoPC(x,y,Z',X)


//A = [[640000,7761.61,70480,1];[1960000,51121.21,316540,1];[4000000,65587.21,512200,1];[6250000,39441.96,496500,1]]
//b = [201.636;188.94;201.588;231.488]

//res = 188+218.5 .*n((0.53 .*n((x-1300)) + 8.0 .*n((y-266).^2)  -2.7 .*n((x-1300).*(y-266))));
