function [Pnew,pnew,rnew,Unew,unew]=Riccatti_contr(M,S,N,m,n,d, ...
                                                   A,B,C,c,Sig, ...
                                                   Theta,theta,P,p,r)

    // -------------------------------------------------------------
    // Recursion de Riccatti sur un pas de temps dans le cas general
    // avec prise en compte d'une contrainte lineaire non dynamique
    // -------------------------------------------------------------

    // Matrices et vecteurs utiles
    E = S' + (B'*P*A);
    F = N + (B'*P*B);
    invF = inv(F);
    G = M + (A'*P*A);
    f = n + (B'*p) + (B'*P*c);
    g = m + (A'*p) + (A'*P*c);
    h = (0.5*trace(Sig*C'*P*C)) + (0.5*c'*P*c) + (p'*c) + d + r;
    invmat = inv(Theta*invF*Theta');
    L = invmat * Theta * invF * E;
    l = invmat * ((Theta*invF*f)+theta);
    K = (Theta'*L) - E;
    k = (Theta'*l) - f;
    // Recursion de Riccatti
    Pnew = G + (2*E'*invF*K) + (K'*invF*K);
    pnew = g + (K'*invF*f) + (K'*invF*k) + (E'*invF*k);
    rnew = h + (0.5*k'*invF*k) + (f'*invF*k);
    // Controle optimal
    Unew = invF * K;
    unew = invF * k;
    //
    // La fonction de Bellman est de la forme
    //     0.5 x Pnew x  +  pnew x  +  rnew
    // et le feedback optimal est de la forme
    //     Unew x  +  unew

endfunction
