function Xn = RacineTchebychev(n,a,b)
    i = linspace(1,n,n)
    Xn=cos(((2.*i-1).*%pi)./(2. *n))
    Xn = (Xn.*((a-b)/2))+((a+b)/2)
endfunction

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
    for i = linspace(0,n-1,n)
            res1(i+1) = 800+((1400/(n+3))*(i+2));
            res2(i+1) = K(i+1);
    end
    res = [res1 res2];
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
p = 600;
//afficheCourbe1(matPuissMdreCre,nbPoints)
ech = linspace(200,2500,p);
couple = ((30/%pi)*fMC(ech,X2))./ech;
puissance = fMC(ech,X2);


//--------------------------------
function res = n(x)
    res = x ./max(x);
endfunction
function res = ff(x,y,A)
    res = (A(1)*(x.^2))+(A(2)*(x.^-1))+(A(3)*(y))+(A(4)*(y.*x))+A(5);
endfunction


function res = gg2(Z,a,M,alpha)
    for i = linspace(1,size(Z,1),size(Z,1))
        for j = linspace(1,size(Z,2),size(Z,2))
            Z(i,j) = Z(i,j)*(1+alpha*(1-((M(1)*(a(i)^2)+M(2)*a(i)+M(3))/Z(i,j))));
        end
    end
    res = Z;
endfunction

function res = matriceVal3D(x,mcP,mcC)
    p = fMC(x,mcP);
    p = p';
    c = fMC(x,mcC);
    c=c';
    res = [(x.^2)' (x.^-1)' p (p.*x') ones(size(x,2),1)];
    res = inv(res'*res)*res'*c;
endfunction

function res = afficheConsoPC(x,y,Z,X)
    res = 0;
    for i = (1:3:p)
        a = fMC(x(i),X);
            q = 0;
        for j = (1:2:p)
            if((Z(i,j) >= a-0.5) & (Z(i,j) <= a+0.5) & q < 5) then
                plot(x(i),y(j),'x');
                res = res + 1;
                q = q+1;
                j = j+1;
            end
        end
    end
endfunction
rpm = linspace(200,2500,p)
c = linspace(1,1700,p);
x = rpm;
y = fMC(x,X2);
a = RacineTchebychev(51,600,2450);
//a = linspace(400,2500,200);

[A,B] = meshgrid(x,y);
Z = ff(A,B,matriceVal3D(a,X2,X));
Z = Z';
//Z = gg2(Z,rpm,X,1.2)
f=gcf();f.color_map=hotcolormap(22);xtitle("interpolation BSFC : f(x,y) = b*x^2 + c*x + f*y + g*y*x + h","regime moteur (tr/min)","puissance fourni (kW)")
zoom_rect([500 0 2600 300]);
colorbar(min(Z),(max(Z)));
grayplot(x,y,Z);
plot(rpm,puissance)
r = afficheConsoPC(x,y,Z,X)


//A = [[640000,7761.61,70480,1];[1960000,51121.21,316540,1];[4000000,65587.21,512200,1];[6250000,39441.96,496500,1]]
//b = [201.636;188.94;201.588;231.488]

//res = 188+218.5 .*n((0.53 .*n((x-1300)) + 8.0 .*n((y-266).^2)  -2.7 .*n((x-1300).*(y-266))));
