// ============================================================
// Solution du probleme global LQG par Riccatti avec contrainte
// ============================================================

clear;

// ------------------
// Donnees generales
// -----------------

// Numero de fenetre graphique
numwin = 0;
// Generateur aleatoire
grand("setsd",123456789);
// Fonctions utiles
exec('Data/Get_data_prob.sce');
exec('Data/Get_data_clas.sce');
exec('Riccatti_contr.sce');
exec('Simulation_bell_temps.sce');
exec('Visualisation.sce');
exec('True_upper_bound.sce');

// ------------------------------
// Donnees du probleme global LQG
// ------------------------------

// Donnees du probleme initial
exec('Data/Data_global.sce');

// Scission des blocs temporels initiaux
exec('Data/Data_global_scission.sce');
// - coefficient de scission des blocs initiaux
//   (valeurs possibles : 1, 2, 3, 4, 6, 8, 12 ou 24)
coefscis = 1;
// - scission proprement dite
[D,TD,Ncl,Acl,TC,TDi,TDf,classes,Data_glob]= ...
      Data_global_scission(D,TD,Ncl,Acl,TC,TDi,TDf,classes,coefscis);

// ---------------------------------------
// Resolution par Riccatti avec contrainte
// ---------------------------------------

// On traite le probleme LQG en resolvant l'equation
// de Riccatti dans laquelle la contrainte couplante
// est prise en compte.

// Initialisation des fonctions de Bellman et des feedbacks optimaux
Bqua = list();
Blin = list();
Bcon = list();
Feed = list();
feed =list();
for t=1:T+1
    Bqua(t) = 0;
    Blin(t) = 0;
    Bcon(t) = 0;
end
for t=1:T
    Feed(t) = 0;
    feed(t) = 0;
end
// Fonction de Bellman a l'instant final
P = MT;
p = mT;
r = dT;
Bqua(T+1) = P;
Blin(T+1) = p;
Bcon(T+1) = r;
// Recursion de Riccatti
for d=D:-1:1
    // - donnees propres au bloc temporel
    icl = Acl(d);
    Data = classes(icl);
    [Mt,Nt,St,mt,nt,dt,At,Bt,Ct,ct,Moy,Sig,Theta,theta]= Get_data_clas(Data);
    // - recursion temporelle
    for h=TD(d):-1:1
        t = TDi(d) + h - 1;
        // . calcul de Riccatti sur un pas de temps
        [Pnew,pnew,rnew,Unew,unew]=Riccatti_contr(Mt,St,Nt, ...
                                                  mt(:,h),nt(:,h),dt(:,h), ...
                                                  At,Bt,Ct,ct(:,h),Sig, ...
                                                  Theta,theta(:,h),P,p,r);
        // . stockage des fonctions de Bellman et des feedbacks
        Bqua(t) = Pnew;
        Blin(t) = pnew;
        Bcon(t) = rnew;
        Feed(t) = Unew;
        feed(t) = unew;
        // . mise a jour pour le pas de temps precedent
        P = Pnew;
        p = pnew;
        r = rnew;
    end
end
// Cout de Bellman a l'instant initial
costBel = (0.5*xini'*P*xini) + (p'*xini) + r;

// -----------------
// Simulation finale
// -----------------

// On effectue la simulation a l'aide des fonctions de Bellman
// obtenues lors de la recursion de Riccatti : on enchaine les
// pas de temps le long de scenarios de bruit, en resolvant un
// probleme ou l'on minimise la somme du cout instantane et du
// cout de Bellman au pas de temps suivant, en propageant d'un
// pas de temps a l'autre l'etat du systeme
//
// Parametre de l'affichage
viewsim = 2;
// Parametres de Monte Carlo
Nsim = 1000001;
Ninc = 50000;
// Simulation proprement dite
grand("setsd",123456789);
[costSim,listG,listg,numwin]=Simulation_bell_temps(numwin,viewsim, ...
                             Nsim,Ninc,costBel,Bqua,Blin,Bcon, ...
                             T,D,TD,Acl,nx,nu,nc,xini,MT,mT,dT,classes);
// Calcul de la borne superieure exacte
[upperBound]=True_upper_bound(listG,listg, ...
                              T,D,TD,Acl,nx,nu,nc,xini,MT,mT,dT,classes);

// ===========
// Fin du code
// ===========
