function [costSim,listG,listg,numwin]=Simulation_bell_temps(numwin,viewsim, ...
                                      Nsim,Ninc,costBel,Bqua,Blin,Bcon, ...
                                      T,D,TD,Acl,nx,nu,nc,xini,MT,mT,dT,classes)

// ----------------------------------------------------
// Simulation du probleme LQG avec contrainte couplante
// ----------------------------------------------------

// Le calcul du feedback optimal sur un pas de temps se fait
// en minimisant la somme du cout instantane et du cout futur
// donne par la fonction de Bellman du pas de temps suivant,
// sous la contrainte couplante

// Initialisation
listu = list();
listx = list();
listw = list();
listG = list();
listg = list();
trajc = zeros(Nsim,T+1);
ctt = [];

// Boucle sur les temps
x = diag(xini) * ones(nx,Nsim);
listx($+1) = x;
t = 0;
for d=1:D
    icl = Acl(d);
    Data = classes(icl);
    [Mt,Nt,St,mt,nt,dt,At,Bt,Ct,ct,Moy,Sig,Theta,theta]= Get_data_clas(Data);
    ctt = [ ctt ; ct' ];
    for h=1:TD(d)
        t = t + 1;
        //  - fonction de Bellman
        P = Bqua(t+1);
        p = Blin(t+1);
        r = Bcon(t+1);
        // - calcul du control optimal
        u = zeros(nu,Nsim);
        E = St' + (Bt'*P*At);
        F = Nt' + (Bt'*P*Bt);
        invF = inv(F);
        f = nt(:,h) + (Bt'*p) + (Bt'*P*ct(:,h));
        invT = inv(Theta*invF*Theta');
        L = invT * Theta * invF * E;
        l = invT * ((Theta*invF*f)+theta(:,h));
        lambda = (L*x) + (diag(l)*ones(nc,Nsim));
        u = - (invF*E*x) - (diag(invF*f)*ones(nu,Nsim)) + (invF*Theta'*lambda);
        // - matrice et vecteur du feedback
        K = (Theta'*L) - E;
        k = (Theta'*l) - f
        G = invF * K;
        g = invF * k;
        listG($+1) = G;
        listg($+1) = g;
        // - bruit sur la dynamique
        w = grand(Nsim,"mn",Moy,Sig);
        // - etat suivant
        xnew = (At*x) + (Bt*u) + (Ct*w) + (diag(ct(:,h))*ones(nx,Nsim));
        // - cout instantane
        trajc(:,t) = sum(0.5*([x;u].*([Mt St;St' Nt]*[x;u])),1) + ...
                     (mt(:,h)'*x) + (nt(:,h)'*u) + dt(:,h);
        // - trajectoires
        listu($+1) = u;
        listx($+1) = xnew;
        listw($+1) = w;
        // - mise a jour de l'etat
        x = xnew;
    end
end
// Cout final
trajc(:,T+1) = sum(0.5*(x.*(MT*x)),1) + (mT'*x) + dT;
// Cout moyen de la simulation
costSim = sum(trajc) / Nsim;
// Ecart-type experimental de ce cout
ecart = sum(trajc,2) - costSim;
ecttSim = sqrt(sum(ecart.*ecart)/(Nsim-1));
// Intervalle de confiance 95%
IClow = costSim - ((1.96*ecttSim)/sqrt(Nsim));
IChigh = costSim + ((1.96*ecttSim)/sqrt(Nsim));
// Affichage des couts
if viewsim>=1 then
    printf('\n Pb LQG - Valeur initial de Bellman : %f' + ...
           '\n          Simulation de Monte Carlo : %f, IC95 : [%f %f]', ...
           costBel,costSim,IClow,IChigh);
end

// Visualisation de la simulation
[numwin]=Visualisation(numwin,viewsim, ...
                       Nsim,Ninc,T,ctt,listx,listu,listw,trajc);

endfunction
