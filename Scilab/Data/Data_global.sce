// -----------------------------------------------------------------------
// Donnees du probleme global avec contrainte couplante et blocs temporels
// -----------------------------------------------------------------------

// Horizon et blocs temporels
// --------------------------

// Horizon d'optimisation
T = 72;
// Nombre de blocs temporels
D = 3;
// Duree des blocs
TD = [ 24 24 24 ];
// Nombre de classes de periodicite
Ncl = 1;
// Affectation des blocs aux classes
Acl = [ 1 1 1 ];
// Duree des blocs d'une meme classe
TC = zeros(1,Ncl);
for icl=1:Ncl
    Tcl = 0;
    for d=1:D
        if Acl(d)==icl then
            if Tcl==0 then
                Tcl = TD(d);
            else
                if Tcl~=TD(d) then
                    disp('Data_Global : durees de bloc incompatibles!');
                    pause
                end
            end
        end
    end
    TC(icl) = Tcl;
end
// Instants initial et final des blocs
TDi = zeros(D,1);
TDf = zeros(D,1);
t = 0;
for d=1:D
    TDi(d) = t +1;
    TDf(d) = t + TD(d);
    t = t + TD(d);
end

// Dimensions du probleme et des sous-problemes
// --------------------------------------------

// Dimensions du probleme global (etat, commande, contrainte)
nx = 4;
nu = 4;
nc = 1;
// Nombre de sous-problemes et dimensions
NLQG = 3;
nxss = [ 1 2 1 ];
nuss = [ 1 2 1 ];
// Etat initial
xini = [ 2.0 ; 3.0 ; 4.0 ; 5.0 ];
// Cout final
coefMT = 1.0;
MT = coefMT * ...
     [ 1.0  0.0  0.0  0.0 ; ...
       0.0  1.0  0.0  0.0 ; ...
       0.0  0.0  1.0  0.0 ; ...
       0.0  0.0  0.0  1.0 ];
mT = coefMT * [ -2.0 ; -3.0 ; -4.0 ; -5.0 ];
dT = [ 0.0 ];

// Classes de periodicite
// ----------------------

classes = list();
// Classe No. 1
icl = 1;
// - cout instantane
coefMt = 2.0;
Mt = coefMt * [ 0.10  0.00  0.00  0.00 ; ...
                0.00  0.10  0.00  0.00 ; ...
                0.00  0.00  0.10  0.00 ; ...
                0.00  0.00  0.00  0.10 ];
coefNt = 2.0;
Nt = coefNt * [ 0.10  0.00  0.00  0.00 ; ...
                0.00  0.10  0.00  0.00 ; ...
                0.00  0.00  0.10  0.00 ; ...
                0.00  0.00  0.00  0.10 ];
St = [ 0.0   0.0   0.0   0.0 ; ...
       0.0   0.0   0.0   0.0 ; ...
       0.0   0.0   0.0   0.0 ; ...
       0.0   0.0   0.0   0.0 ];
mt = [ 0.0 ; 0.0 ; 0.0 ; 0.0 ];
mt = diag(mt) * ones(nx,TC(icl));
nt = [ -1.0 ; -1.0 ; -1.0 ; -1.5 ];
nt = diag(nt) * ones(nu,TC(icl));
nt = [ -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 ... // u1
       -0.90 -0.80 -0.70 -0.60 -0.60 -0.70 -0.80 -0.90 ...
       -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 ;
       -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 ... // u2
       -0.95 -0.90 -0.85 -0.80 -0.80 -0.85 -0.90 -0.95 ...
       -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 ;
       -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 ... // u3
       -0.95 -0.90 -0.85 -0.80 -0.80 -0.85 -0.90 -0.95 ...
       -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 ;
       -1.50 -1.50 -1.50 -1.50 -1.50 -1.50 -1.50 -1.50 ... // u4
       -1.55 -1.60 -1.65 -1.70 -1.70 -1.65 -1.60 -1.55 ...
       -1.50 -1.50 -1.50 -1.50 -1.50 -1.50 -1.50 -1.50 ];
// Terme correspondant à la moyenne de la demande
// ajuste pour que le cout soit de l'ordre de 100
dt = [ 2.2005 ];
dt = diag(dt) * ones(1,TC(icl));
// - dynamique
At = [ 1.0  0.0  0.0  0.0 ; ...
       0.0  1.0  0.0  0.0 ; ...
       0.0  0.0  1.0  0.0 ; ...
       0.0  0.0  0.0  1.0 ];
Bt = [ -1.0   0.0   0.0   0.0 ; ...
        0.0  -1.0   0.0   0.0 ; ...
        0.0   1.0  -1.0   0.0 ; ...
        0.0   0.0   0.0  -1.0 ];
Ct = [ 1.0   0.0   0.0   0.0 ; ...
       0.0   1.0   0.0   0.0 ; ...
       0.0   0.0   1.0   0.0 ; ...
       0.0   0.0   0.0   1.0 ];
ct = [ 1.0 ; 1.0 ; 1.0 ; 1.0 ];
ct = diag(ct) * ones(nx,TC(icl));
//  - bruit
Moy = [ 0.0 ; 0.0 ; 0.0 ; 0.0 ];
Sig = [ 0.01  0.0   0.0  0.0 ; ...
        0.0   0.01  0.0  0.0 ; ...
        0.0   0.0  0.01  0.0 ; ...
        0.0   0.0   0.0  0.01 ];
// - contrainte
Theta = [ 1.0 1.0 1.0 1.0 ];
theta = [ 5.0 ];
theta = diag(theta) * ones(nc,TC(icl));
// - ensemble des donnees de la classe
classes($+1) = list(Mt,Nt,St,mt,nt,dt,At,Bt,Ct,ct,Moy,Sig,Theta,theta);

// Stockage des donnees du probleme global
Data_glob = list(nx,nu,nc,xini,MT,mT,dT,classes);

