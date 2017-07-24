//Calcul la solution au sens des moindres carrés
//x : un vecteur de points (données)
//y : un vecteur de points (valeur) ax^n+bx^(n-1)+..p = y
//n : l'ordre du modèle
function X = moindresCarres(x,val,ordre)
    A = [x,ones(size(x,1),1)];
    xT = x;
    for i = (2:1:ordre)
        xT = xT.*x;
        A = [xT A];
    end
    X = inv(A'*A)*A'*val;
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

//création une matrice [x y] de n éléments ; pour le calcul des moindres carrés
//démarre à 1000, echantillonage tout les 100. 
function res = donnees(n,K)
    for i = linspace(0,n-1,n)
            res1(i+1) = 800+((1400/(n+3))*(i+2));
            res2(i+1) = K(i+1);
    end
    res = [res1 res2];
endfunction


//Acquisition des points
nbPoints = 11;
nbPtsConso  = 10;
//Points pour moindres carrés
ptsPuiss = [178,196,214,230,246,261,263,264,265,264,261];
p1 = [193,190,189,188,189,191,193,195,198,201];
p2 = [150,208,203,195,155,197,200,206,210,202];
p3 = [210,213,209,201,206,207,210,217,225,239];
p4 = [240,234,229,215,221,230,237,240,250,260];
p5 = [260,256,249,240,242,217,250,259,270,285];
p6 = [280,278,269,266,269,275,282,290,300,310];
p7 = [300,297,289,280,283,290,299,310,320,330];
p8 = [340,336,320,315,317,320,332,345,350,380];

matConso = [p1;p2;p3;p4;p5];
RPM = 0;
for i = linspace(0,nbPtsConso-1,nbPtsConso)
    RPM(i+1) = 800+((1400/(nbPtsConso+3))*(i+2));
end


matPuissMdreCre = donnees(nbPoints,ptsPuiss);

x2 = matPuissMdreCre(:,1);
y2 = matPuissMdreCre(:,2);
degre = 2;
X = [];
for i = linspace(1,5,5)
    X = [X;moindresCarres(RPM,(matConso(i,:))',degre)'];
end
X2 = moindresCarres(x2,y2,degre);

p = 600;
//afficheCourbe1(matPuissMdreCre,nbPoints)
ech = linspace(50,2500,p);

couple = ((30/%pi)*fMC(ech,X2))./ech;
puissance = fMC(ech,X2);


//--------------------------------

function res = calculGrilleSurface(x,y,A)
    res = 0;
    for i = linspace(0,3,4)
        res = res + A(2*i+1)*(x.^(4-i)) + A(2*(i+1))*(y.^(4-i));
    end
    res = res + A(9)*((x).*(y));
    res = res + A(10) + A(11)*((x+1).^-1);
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
    x = []
    for i = linspace(1,t,t)
        x = [x a];
    end
    fMCa = fMC(a,mcP);
    p = [fMCa];
    for i = linspace(t-1,1,t-1)
        p = [p (i/t)*fMCa];
    end
    p = p';
    c = [fMC(a,MC(1,:)')];
    for i = linspace(2,t,t-1)
        c = [c fMC(a,MC(i,:)')];
    end
    c=c';
    res = [x'.^4 p.^4];
    for i = linspace(1,3,3)
        res = [res x'.^(4-i) p.^(4-i)];
    end
    a = (x').*(p);
    res = [res a];
    res = [res ones(size(x,2),1) ];
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
    for i = (1:7:p)
        a = fMC(x(i),X);
            q = 0;
        for j = (1:7:p)
            if((Z(i,j) >= a-3.2) & (Z(i,j) <= a+3.2) & q < 20) then
                plot(x(i),y(j),'x');
                res = res + 1;
                q = q+1;
            end
        end
    end
endfunction
rpm = linspace(50,2500,p);
c = linspace(1,1700,p);
x = rpm;
y = fMC(x,X2);
a = linspace(400,2500,50);

[A,B] = meshgrid(x,y);
Z = calculGrilleSurface(A,B,matriceVal3D(size(X,1),a,X2,X));
Z = Z';
//Z = gg2(Z,rpm,X,1.2)
f=gcf();f.color_map=hotcolormap(36);xtitle("Graphique d interpolation d un BSFC diesel : f(x,y)=ax²+by+c","regime moteur (tr/min)","puissance fourni (kW)")
zoom_rect([500 0 2600 300]);
colorbar(150,380);
grayplot(x,y,Z);
plot(rpm,puissance)
r = afficheConsoPC(x,y,Z,X(1,:)')
r = afficheConsoPC(x,y,Z,X(3,:)')
