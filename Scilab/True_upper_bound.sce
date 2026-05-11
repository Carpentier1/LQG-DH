function [upperBound]=True_upper_bound(listG,listg, ...
                                       T,D,TD,Acl,nx,nu,nc,xini,MT,mT,dT, ...
                                       classes)

// ----------------------------------------------------
// Calcul d'une borne superieure exacte du cout optimal
// ----------------------------------------------------

// Le calcul se fait par l'equation de programmation dynamique
// sur le probleme global, equation dans  laquelle la commande
// utilisee est celle obtenue lors de la simulation  finale du
// probleme. Cette commande est telle qu'elle verifie toujours
// la contrainte couplante. La simulation finale fournit quant
// a elle une borne superieure statistique du cout optimal.

// Initialisation des fonctions de Bellman
Bqua = list();
Blin = list();
Bcon = list();
for t=1:T+1
    Bqua(t) = 0;
    Blin(t) = 0;
    Bcon(t) = 0;
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
        // . coefficients du controle
        Gt = listG(t);
        gt = listg(t);
        // . modification du cout et de la dynamique
        Mtb = Mt + ((St*Gt)+(Gt'*St')) + (Gt'*Nt*Gt);
        mtb = (St*gt) + (Gt'*Nt*gt) + mt(:,h) + (Gt'*nt(:,h));
        dtb = (0.5*gt'*Nt*gt) + (nt(:,h)'*gt) + dt(:,h);
        Atb = At + (Bt*Gt);
        ctb = ct(:,h) + (Bt*gt);
        // . calcul de Riccatti sur un pas de temps
        Pnew = Mtb + (Atb'*P*Atb);
        pnew = mtb + (Atb'*(p+(P*ctb)));
        rnew = dtb + (0.5*ctb'*P*ctb) + (0.5*trace(Sig*Ct'*P*Ct)) + ...
               (p'*ctb) + r;
        // . stockage des fonctions de Bellman et des feedbacks
        Bqua(t) = Pnew;
        Blin(t) = pnew;
        Bcon(t) = rnew;
        // . mise a jour pour le pas de temps precedent
        P = Pnew;
        p = pnew;
        r = rnew;
    end
end
// Borne superieure exacte
upperBound = (0.5*xini'*P*xini) + (p'*xini) + r;
printf('\n          Borne superieure exacte   : %f',upperBound);

endfunction

