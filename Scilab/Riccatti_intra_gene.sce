function [Pnew,Anew,bnew,Dnew,dnew,enew,Unew,unew,vnew]= ...
               Riccatti_intra_gene(M,S,N,m,n,d,A,B,C,c,Sig, ...
                                   P,Alpha,betaa,Delta,delta,epsilon)

    // -------------------------------------------------------------
    // Recursion de Riccatti sur un pas de temps dans le cas general
    // en l'absence de contrainte lineaire non dynamique, et pour la
    // quelle la dependance lineaire en le multiplicateur nu associe
    // au couplage temporel est calculee de maniere explicite
    // -------------------------------------------------------------

    // Matrices et vecteurs utiles
    E = S' + (B'*P*A);
    F = N + (B'*P*B);
    invF = inv(F);
    G = M + (A'*P*A);
    Alphah = B' * Alpha;
    betaah = n + (B'*betaa) + (B'*P*c);
    Alphat =  A' * Alpha;
    betaat = m + (A'*betaa) + (A'*P*c);
    deltah = delta + (Alpha'*c);
    epsilonh = epsilon + d + (c'*betaa) + (0.5*c'*P*c) + ...
               (0.5*trace(Sig*C'*P*C)); 
    //
    // Lien avec les vecteurs utilises dans la recursion
    // f = Alphah*nu + betaah
    // g = Alphat*nu + betaat
    // h = 0.5*nu*Delta*nu + deltah*nu + epsilonh
    //
    // Recursion de Riccatti
    Pnew = G - (E'*invF*E);
    Anew = Alphat - (E'*invF*Alphah);
    bnew = betaat - (E'*invF*betaah);
    Dnew = Delta - (Alphah'*invF*Alphah);
    dnew = deltah - (Alphah'*invF*betaah);
    enew = epsilonh - (0.5*betaah'*invF*betaah);
    // Controle optimal
    Unew = - invF*E;
    unew = - invF*Alphah;
    vnew = - invF*betaah;
    //
    // La fonction de Bellman est de la forme
    //     0.5 x Pnew x  +  pnew x  +  rnew
    // avec
    //     pnew = Anew nu  +  bnew
    //     rnew = 0.5 nu Dnew nu  +  dnew nu  +  enew
    // et le feedback est de la forme
    //     Unew x  +  unew nu  +  vnew

endfunction
