function [Pnew,Anew,bnew,Dnew,dnew,enew,Unew,unew,vnew]= ...
               Riccatti_intra_contr(M,S,N,m,n,d, ...
                                    A,B,C,c,Sig,Theta,theta, ...
                                    P,Alpha,betaa,Delta,delta,epsilon)

    // -------------------------------------------------------------
    // Recursion de Riccatti sur un pas de temps dans le cas general
    // avec prise en compte d'une contrainte lineaire non dynamique,
    // et ou la dependance lineaire en le multiplicateur nu associe
    // au couplage temporel est calculee de maniere explicite
    // -------------------------------------------------------------

    // Matrices et vecteurs utiles
    E = S' + (B'*P*A);
    F = N + (B'*P*B);
    invF = inv(F);
    G = M + (A'*P*A);
    invmat = inv(Theta*invF*Theta');
    L = invmat * (Theta*invF*E);
    K = (Theta'*L) - E;
    Alphah = B' * Alpha;
    betaah = n + (B'*betaa) + (B'*P*c);
    Alphat =  A' * Alpha;
    betaat = m + (A'*betaa) + (A'*P*c);
    Alphac = invmat * Theta * invF * Alphah;
    betaac = (invmat*Theta*invF*betaah) + (invmat*theta);
    Alphap = (Theta'*Alphac) - Alphah;
    betaap = (Theta'*betaac) - betaah;
    deltah = delta + (Alpha'*c);
    epsilonh = epsilon + d + (c'*betaa) + (0.5*c'*P*c) + ...
               (0.5*trace(Sig*C'*P*C));
    //
    // Lien avec les vecteurs utilises dans la recursion
    // f = Alphah*nu + betaah
    // g = Alphat*nu + betaat
    // l = Alphac*nu + betaac
    // k = Alphap*nu + betaap
    // h = 0.5*nu*Delta*nu + deltah*nu + epsilonh
    //
    // Recursion de Riccatti
    Pnew =G + (2*E'*invF*K) + (K'*invF*K);
    Anew = Alphat + (K'*invF*Alphah) + (K'*invF*Alphap) + (E'*invF*Alphap);
    bnew = betaat + (K'*invF*betaah) + (K'*invF*betaap) + (E'*invF*betaap);
    Dnew = Delta + (Alphap'*invF*Alphap) + (2*Alphah'*invF*Alphap);
    dnew = deltah + (Alphap'*invF*(betaah+betaap)) + (Alphah'*invF*betaap);
    enew = epsilonh + (0.5*betaap'*invF*betaap) + (betaah'*invF*betaap);
    // Controle optimal
    Unew = invF * K;
    unew = invF * Alphap;
    vnew = invF * betaap;
    //
    // La fonction de Bellman est de la forme
    //     0.5 x Pnew x  +  pnew x  +  rnew
    // avec
    //     pnew = Anew nu  +  bnew
    //     rnew = 0.5 nu Dnew nu  +  dnew nu  +  enew
    // et le feedback est de la forme
    //     Unew x  +  unew nu  +  vnew

endfunction
