//Acquisition des points
nbPoints = 11;
nbPtsConso  = 10;
//Points pour moindres carrés
ptsPuiss = [178,196,214,230,246,261,263,264,265,264,261];
ptsConso = [193,190,189,188,189,191,193,195,198,201];


ptsConso1 = [217.3,210.7,205,203,190,189,188,189,191,193,195,198,201,207,215.2];
MatriceConso = []
for i = linspace(1,15,10)
    MatriceConso = [MatriceConso ; ptsConso1];
end
Mtranslation = [[1.027,1.011,1.,1.04,1.01,1,0.99,1,1,1,1,1,1,1,1];[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]
ptsConsoDemi = ptsConso*1.46; //x2
ptsConsoTQ = ptsConso*1.23; //x1.5

function Mres = matriceTranslation(rpmMin,rpmMax,echX,pMin,pMax,echY,MatConso)
    
endfunction
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

//attend 15 points en théorie
function res = d2(n,K)
    for i = linspace(1,n,n)
            res1(i) = 700+(100*i);
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
matConsoUMdreCre = donnees(nbPtsConso,ptsConsoUQ);
matConsoDMdreCre = donnees(nbPtsConso,ptsConsoDemi);
matConsoTMdreCre = donnees(nbPtsConso,ptsConsoTQ);
matPuissMdreCre = donnees(nbPoints,ptsPuiss);

x = matConsoMdreCre(:,1);
xu = matConsoUMdreCre(:,1);
xd = matConsoDMdreCre(:,1);
xt = matConsoTMdreCre(:,1);
x2 = matPuissMdreCre(:,1);
yu = matConsoUMdreCre(:,2);
yd = matConsoDMdreCre(:,2);
yt = matConsoTMdreCre(:,2);
y = matConsoMdreCre(:,2);
y2 = matPuissMdreCre(:,2);
degre = 2;
X = moindresCarres(x,y,degre);
XU = moindresCarres(xu,yu,degre);
XD = moindresCarres(xd,yd,degre);
XT = moindresCarres(xt,yt,degre);
X2 = moindresCarres(x2,y2,degre);

p = 600;
//afficheCourbe1(matPuissMdreCre,nbPoints)
ech = linspace(500,2500,p);

couple = ((30/%pi)*fMC(ech,X2))./ech;
puissance = fMC(ech,X2);


//--------------------------------

function res = calculGrilleSurface(x,y,A)
    res = 0;
    for i = linspace(0,3,4)
        res = res + A(2*i+1)*(x.^(4-i)) + A(2*(i+1))*(y.^(4-i));
    end
    res = res + A(9)*((x).*(y));
    res = res + A(10);
endfunction

function res = gg2(Z,a,M,alpha)
    for i = linspace(1,size(Z,1),size(Z,1))
        for j = linspace(1,size(Z,2),size(Z,2))
            Z(i,j) = Z(i,j)*(1+alpha*(1-((M(1)*(a(i)^2)+M(2)*a(i)+M(3))/Z(i,j))));
        end
    end
    res = Z;
endfunction

function res = matriceVal3D(t,a,mcP,MC)
    x = [a a a a]
    p = [((mcP(1)*(a.^2))+(mcP(2)*a)+mcP(3))]
    for i = linspace(t-1,1,t-1)
        p = [p (i/t)*((mcP(1)*(a.^2))+(mcP(2)*a)+mcP(3))];
    end
    p = p';
    c = [((MC(1,1)*(a.^2))+(MC(1,2)*a)+MC(1,3))]
    for i = linspace(2,t,t-1)
        c = [c (i/t)*((MC(i,1)*(a.^2))+(MC(i,2)*a)+MC(i,3))];
    end
    c=c';
    res = [x'.^4 p.^4];
    for i = linspace(1,3,3)
        res = [res x'.^(4-i) p.^(4-i)];
    end
    a = (x').*(p);
    res = [res a];
    res = [res ones(size(x,2),1)];
    res = inv(res'*res)*res'*c;
endfunction

function res = matriceEch(taille,echX,matricePtsConso)
    res = [];
    for i = (1:1:size(matricePtsConso,1))
        matConso = d2(size(matricePtsConso,2),matricePtsConso(i,:));
        MT = moindresCarres(matConso(:,1),matConso(:,2),2)
        res = [res ; MT];
    end
endfunction
function res = afficheConsoPC(x,y,Z,X)
    res = 0;
    for i = (1:5:p)
        a = fMC(x(i),X);
            q = 0;
        for j = (1:6:p)
            if((Z(i,j) >= a-3.2) & (Z(i,j) <= a+3.2) & q < 20) then
                plot(x(i),y(j),'x');
                res = res + 1;
                q = q+1;
                j = j+2
            end
        end
    end
endfunction
rpm = linspace(500,2500,p)
c = linspace(1,1700,p);
x = rpm;
y = (X2(1)*x.^2+X2(2)*x+X2(3));
a = linspace(400,2500,250);

[A,B] = meshgrid(x,y);
MX = matriceEch(10,a,Points)
Z = calculGrilleSurface(A,B,matriceVal3D(10,a,X2,MX));
Z = Z';
//Z = gg2(Z,rpm,X,1.2)
f=gcf();f.color_map=hotcolormap(32);xtitle("Graphique d interpolation d un BSFC diesel : f(x,y)=ax²+by+c","regime moteur (tr/min)","puissance fourni (kW)")
zoom_rect([500 0 2600 300]);
colorbar(min(Z),(max(Z)));
grayplot(x,y,Z);
plot(rpm,puissance)
r = afficheConsoPC(x,y,Z,X)


//A = [[640000,7761.61,70480,1];[1960000,51121.21,316540,1];[4000000,65587.21,512200,1];[6250000,39441.96,496500,1]]
//b = [201.636;188.94;201.588;231.488]

//res = 188+218.5 .*n((0.53 .*n((x-1300)) + 8.0 .*n((y-266).^2)  -2.7 .*n((x-1300).*(y-266))));
