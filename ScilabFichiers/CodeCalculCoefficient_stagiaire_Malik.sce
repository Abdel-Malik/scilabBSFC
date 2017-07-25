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
ptsCouple = [1696,1700,1696,1690,1680,1659,1578,1493,1415,1334,1253];
ptsConso = [193,190,189,188,189,191,193,195,198,201];
//ptsConsoModifiee = [193,190,189,188,189,191,193,195,198,206];

//exposant pour le calcul de consommation spécifique
expRotation1 = 1.5;
expRotation2 = -1;
expCorrelation = 0.9;
expCouple = 0.45;

/*Fonctions maths*/
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
/*fin fonction*/

//Vecteurs : plage de rotations moteurs en fonction d'un intervalle donné
rpmP = linspace(intervalleBasP,intervalleHautP,size(ptsPuiss,2));
rpmCouple = linspace(intervalleBasC,intervalleHautC,size(ptsCouple,2));
rpmConso = linspace(intervalleBasConso,intervalleHautConso,size(ptsConso,2));
rpmConsoM = linspace(intervalleBasConso,intervalleHautConso,size(ptsConsoModifiee,2));

//--Calcul de coefficients par la méthode des moindres carrés--//
polyCouple = moindresCarres(rpmCouple',ptsCouple',degreInterpolationCouple);
polyConso = moindresCarres(rpmConso',ptsConso',degreInterpolationConso);
polyConsoM = moindresCarres(rpmConsoM',ptsConsoModifiee',degreInterpolationConso);


//vecteur (abscisses) pour le calcul et l'affichage des courbes
ech=linspace(miniR,maxiR,1000);

//Affichage Courbes de couple
y3 = fMC(ech,polyCouple);
plot(ech,y3,'c');
plot(rpmCouple,ptsCouple,'b--');
xgrid(1);
zoom_rect([miniR miniCouple maxiR 1.1*maxiCouple]);
xtitle("Courbe de Couple pleine charge","regime moteur (tr/min)","couple moteur (Nm)");
//attends un clique souris pour continuer le code
xclick();
//remise à zéro de l'affichage
clf();

//Affichage Courbes de puissance
plot(ech,(%pi/30)*(y3.*ech),'g');
plot(rpmP,ptsPuiss,'b--');
xtitle("Courbe de puissance pleine charge","regime moteur (tr/min)","puissance (W)");
zoom_rect([miniR miniP maxiR 1.1*maxiP]);
xgrid(1);
xclick();
clf();

//Affichage Courbes de consommation
plot(ech,fMC(ech,polyConso),'r');
plot(rpmConso,ptsConso,'b--');
xgrid(1);
xtitle("Courbe de consommation pleine charge","regime moteur (tr/min)","Consommation (g/kWh)");
zoom_rect([miniR miniConso maxiR 1.1*maxiConso]);
xclick();
clf();

//--Partie graphique consommation spécifique--//

function res = calculGrilleSurface(x,y,A)
    res = A(1)*(x.^expRotation1)+A(2)*((x+1).^expRotation2) + A(3)*(y.*x).^expCorrelation + A(4)*(y.^expCouple) + A(5);
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
    res = [(xx.^expRotation1)' ((xx+1).^expRotation2)' ((t.*xx).^expCorrelation)' (t.^expCouple)' ones(size(x,2),1)];
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

//nombre de points par axe pour le graphique (influe sur le niveau de détail)
ptsGraph = 150;

//initialisation vecteurs
rpm = linspace(400,maxiR,ptsGraph);
couple = fMC(rpm,polyCouple);
puissance =(%pi/30)*(couple.*rpm);
a = linspace(miniR,maxiR,n);
coupleVal = linspace(miniCouple,1.1*max(couple),ptsGraph);
puissVal = linspace(miniP,1.1*max(puissance),ptsGraph);

[A,B] = meshgrid(rpm,coupleVal);

Z = calculGrilleSurface(A,B,matriceVal3D(a,polyCouple,polyConso));
Z = Z';
//Z = amplifierEcart(Z,rpm,polyConso,1.10);

//--Affichage--//

//positionnement de l'affichage
zoom_rect([miniR miniCouple maxiR 1.1*maxiCouple]);

//préparation de la coloration (choix d'un nombre de nuances)
f=gcf();f.color_map=hotcolormap(32);

//nomme le graphique ainsi que les axes
xtitle("Graphique d interpolation d un BSFC diesel : f(x,y)=ax^"+string(expRotation1)+" + bx^("+string(expRotation2)+") + cxy^"+string(expCorrelation)+" + dy^"+string(expCouple)+" + e","regime moteur (tr/min)","couple fourni (Nm)")
//initialise les extrémités de l'intervalle de la matrice utilisé pour la coloration
colorbar(min(Z),(max(Z)));

//coloration (affichage du 3e axe)
grayplot(rpm,coupleVal,Z);

//affichage courbes
plot(rpm,couple);

//affichage points de donnée
plot(rpmCouple,ptsCouple,'roo');

//affichage des croix bleus
r = afficheConsoPC(rpm,coupleVal,Z,polyConso);

//Ecriture coefficients dans la console
disp(polyConso,"Les coefficients du polynôme de consommation (du degré n au degré 0) :");

disp(polyCouple,"Les coefficients du polynôme de couple (du degré n au degré 0):");

disp(matriceVal3D(a,polyCouple,polyConso),"Les coefficients de la formule de consommation spécifique :");
