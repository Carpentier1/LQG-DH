function [numwin]=Visualisation(numwin,viewsim, ...
                                Nsim,Ninc,T,ctt,listx,listu,listw,trajc)

// --------------------------------------------
// Visualisation des resultats d'une simulation
// --------------------------------------------

if viewsim>=2 then
    // Trajectoires des etats
    numwin = numwin + 1;
    scf(numwin);
    clf();
    for nsim=1:Ninc:Nsim
        trajx = [];
        for t=1:T+1
            xt = listx(t);
            trajx = [ trajx ; xt(:,nsim)' ];
        end
        plot2d([0:T],trajx);
    end
    xgrid;
    xtitle('Etats');
    // Trajectoires des commandes
    numwin = numwin + 1;
    scf(numwin);
    clf();
    for nsim = 1:Ninc:Nsim
        traju = [];
        for t=1:T
            ut = listu(t);
            traju = [ traju ; ut(:,nsim)' ];
        end
        plot2d([0:T-1],traju);
    end
    xgrid;
    xtitle('Commandes');
    // Trajectoires des bruits
    numwin = numwin + 1;
    scf(numwin);
    clf();
    for nsim = 1:Ninc:Nsim
        trajw = [];
        for t=1:T
            wt = listw(t);
            trajw = [ trajw ; wt(:,nsim)' ];
        end
        plot2d([1:T],ctt+trajw);
    end
    xgrid;
    xtitle('Aleas');
    // Trajectoires des couts
    numwin = numwin + 1;
    scf(numwin);
    clf();
    for nsim = 1:Ninc:Nsim
        plot2d([0:T],trajc(nsim,:));
    end
    xgrid;
    xtitle('Couts');
end

endfunction
