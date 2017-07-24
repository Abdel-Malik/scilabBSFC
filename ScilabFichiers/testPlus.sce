//**Données du moteur**//
//Plage de fonctionnement (tr/min): 
miniR = 800;
maxiR = 2500;
//Couple fourni (Nm): 
miniCouple = 0;
maxiCouple = 1800;
//Puissance fourni (W): 
miniP = 0;
maxiP = 266000;
//Consommation (g/kWh): 
miniConso = 180;
maxiConso = 210;
n = 10;
//intervalle dans lequel l'échantillonnage a été réalisé
//Puissance :
intervalleBasP = 1000;
intervalleHautP = 2000;
//Couple :
intervalleBasC = 1000;
intervalleHautC = 2000;
//Conso :
intervalleBasConso = 1015.3846;
intervalleHautConso = 1984.6154;

//interpolation de l'echantillonnage
degreInterpolationCouple = 3;
degreInterpolationConso = 2;
//Points échantillonnés (pas constant) pour moindres carrés
ptsPuiss = [172000,194200,212100,229000,246000,261900,263800,264600,265000,264600,260200];
p1 = [1696,1700,1696,1690,1680,1659,1578,1493,1415,1334,1253];
p2 = [193,190,189,188,189,191,193,195,198,201];
p2 = [150,208,203,195,155,197,200,206,210,202];
p3 = [210,213,209,201,206,207,210,217,225,239];
p4 = [240,234,229,215,221,230,237,240,250,260];
p5 = [260,256,249,240,242,217,250,259,270,285];
p6 = [280,278,269,266,269,275,282,290,300,310];
p7 = [300,297,289,280,283,290,299,310,320,330];
p8 = [340,336,320,315,317,320,332,345,350,380];
function Xn = RacineTchebychev(n,a,b)
    i = linspace(1,n,n)
    Xn=cos(((2.*i-1).*%pi)./(2. *n))
    Xn = (Xn.*((a-b)/2))+((a+b)/2)
endfunction

//*Résolution du système des moindres carrés
//retourne X (un vecteur) dans l'équation : AX = b
function X = moindresCarres(x,val,ordre)
    A = [x,ones(size(x,1),1)];
    xT = x;
    for i = (2:1:ordre)
        xT = xT.*x;
        A = [xT A]
    end
    X = inv(A'*A)*A'*val;
endfunction

//x un vecteur de points
//X un vecteur contenant des coefficients des moindres carrés
//y le calcul du polynôme de coefficient X pour tout élément de x
function y = fMC(x,X)
    l = size(X,1)
    y = zeros(1,size(x,2))
    xT = ones(1,size(x,2))
    for i = linspace(l,1,l)
        y = y+xT*X(i);
        xT = xT.*x;
    end
endfunction

//fMC pour Matrice
function y = fMCM(x,X)
    l = size(X,1);
    lM = size(x,1);
    y = [];
    yt = zeros(1,size(x,2));
    xT = ones(1,size(x,2));
    for j = linspace(1,lM,lM)
        for i = linspace(l,1,l)
            yt = yt+xT*X(i);
            xT = xT.*x(j,:);
        end
        y = [y;yt];
    end
endfunction
function res = calculRxy(x,y)
    nb = size(x,2);
    xy = x*y';
    xySep = (sum(x)*sum(y))/nb;
    varX = sum(x.^2)-((sum(x)^2)/nb);
    varY = sum(y.^2)-((sum(y)^2)/nb);
    res = (xy-xySep)/(varX*varY);    
endfunction

//calcul la régression [0;1] entre x et y
function res = calculRcarre(x,y)
    nb = size(x,2);
    xm = sum(x)/nb;
    ym = sum(y)/nb
    Sxy = (sum((x-xm).*(y-ym)))/(nb-1);
    Sx2 = sum((x-xm).^2)/(nb-1);
    Sy2 = sum((y-ym).^2)/(nb-1);
    res = (Sxy^2)/(Sx2*Sy2);    
endfunction

rpmP = linspace(intervalleBasP,intervalleHautP,size(ptsPuiss,2));
rpmCouple = linspace(intervalleBasC,intervalleHautC,size(ptsCouple,2));
rpmConso = linspace(intervalleBasConso,intervalleHautConso,size(ptsConso,2));
polyCouple = moindresCarres(rpmCouple',ptsCouple',degreInterpolationCouple);

//--Etape 1 Calcul de coefficients par la méthode des moindres carrés--//
matConso = [p1;p2;p3;p4;p5];
RPM = 0;
for i = linspace(0,nbPtsConso-1,nbPtsConso)
    RPM(i+1) = 800+((1400/(nbPtsConso+3))*(i+2));
end
X = [];
for i = linspace(1,5,5)
    X = [X;moindresCarres(RPM,(matConso(i,:))',degre)'];
end
ptsGraph = 300;
rpm = linspace(200,maxiR,ptsGraph);
couple = fMC(rpm,polyCouple);
puissance =(%pi/30)*(couple.*rpm);


function res = calculGrilleSurface(x,y,A)
    res = A(1)*(x.^3)+A(2)*((x+1).^-1) + A(3)*(y.*x) + A(4)*(y.^1) + A(5);
endfunction

function res = amplifierEcart(Z,a,M,alpha)
    valeurMin = fMC(a,M);
    consoT = fMC(a,M);
    for regime = linspace(1,size(Z,1),size(Z,1))
        for j = linspace(1,size(Z,2),size(Z,2));
            Z(regime,j) = Z(regime,j)*(1+(alpha*(1-(valeurMin(regime)/Z(regime,j)))));
        end
    end
    res = Z;
endfunction

function res = matriceVal3D(x,mcCouple,mcConso,mMinConso)
    res = [];
    p=0;
    xx =  x;
    t = fMC(x,mcCouple);
    c = fMC(x,mcConso);
    c = c';
    res = [(xx.^3)' ((xx+1)').^-1 (t.*xx)' (t.^1)' ones(size(x,2),1)];
    res = inv(res'*res)*res'*c;
endfunction

function res = afficheConsoPC(x,y,Z,X)
    res = 0;
        XIntervalleConsoMin = [];
        YIntervalleConsoMin = [];
        xIntervalleValeursTh = [];
        yIntervalleValeursTh = [];
    for i = (1:1:ptsGraph)
        a = fMC(x(i),X);
        q = 0;
        for j = (1:2:ptsGraph)
            if(modulo(i,3) == 1) then
                if((Z(i,j) >= a-0.3) & (Z(i,j) <= a+0.3) & q < 20) then
                    XIntervalleConsoMin = [XIntervalleConsoMin x(i)];
                    YIntervalleConsoMin = [YIntervalleConsoMin y(j)];
                    //q = q+1; //A décommenter si le nombre de croix est trop important
                end
            end
            for k = (1:1:size(rpmConso,2))
                p=0;
                if(x(i) <= (rpmConso(k)+5) & x(i) >= (rpmConso(k)-5)) then
                    for jj = (1:1:ptsGraph)
                        if((Z(i,j) >= ptsConso(k)-0.2) & (Z(i,j) <= ptsConso(k)+0.2) & p == 0) then
                            xIntervalleValeursTh = [xIntervalleValeursTh x(i)];
                            yIntervalleValeursTh = [yIntervalleValeursTh y(j)];
                            p = 1;
                        end
                    end
                end
            end
        end
    end
    if(size(xIntervalleValeursTh,2) > 0 & size(yIntervalleValeursTh,2) > 0) then
        plot(xIntervalleValeursTh,yIntervalleValeursTh,'wxx');
    end
    if(size(XIntervalleConsoMin,2) > 0 & size(YIntervalleConsoMin,2) > 0) then
        plot(XIntervalleConsoMin,YIntervalleConsoMin,'x');
    end
endfunction


a = linspace(miniR,maxiR,n);
coupleVal = linspace(miniCouple,1.1*max(couple),ptsGraph);
puissVal = linspace(miniP,1.1*max(puissance),ptsGraph);

[A,B] = meshgrid(rpm,coupleVal);

Z = calculGrilleSurface(A,B,matriceVal3D(a,polyCouple,polyConso));
Z = Z';
//Z = amplifierEcart(Z,rpm,polyConso,1.10);

//Affichage

zoom_rect([miniR miniCouple maxiR 1.1*maxiCouple]);
f=gcf();f.color_map=hotcolormap(36);
xtitle("Graphique d interpolation d un BSFC diesel : f(x,y)=ax^3 + bx^(-1) + cxy + d","regime moteur (tr/min)","couple fourni (Nm)")
colorbar(min(Z),(max(Z)));
grayplot(rpm,coupleVal,Z);
plot(rpm,couple);
plot(rpmCouple,ptsCouple,'roo');
r = afficheConsoPC(rpm,coupleVal,Z,polyConso);

//Ecriture coefficients dans la console
disp(polyConso,"Les coefficients du polynôme de consommation (du degré n au degré 0) :");

disp(polyCouple,"Les coefficients du polynôme de couple (du degré n au degré 0):");

disp(matriceVal3D(a,polyCouple,polyConso),"Les coefficients de la formule de consommation spécifique :");
