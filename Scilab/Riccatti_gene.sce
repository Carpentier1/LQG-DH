function [Pnew,pnew,rnew,Unew,unew]= ...
               Riccatti_gene(M,S,N,m,n,d, ...
                             A,B,C,c,Sig,P,p,r)

    // -------------------------------------------------------------
    // Recursion de Riccatti sur un pas de temps dans le cas general
    // en l'absence de contrainte lineaire non dynamique
    // -------------------------------------------------------------

    // Matrices et vecteurs utiles
    E = S' + (B'*P*A);
    F = N + (B'*P*B);
    invF = inv(F);
    G = M + (A'*P*A);
    f = n + (B'*p) + (B'*P*c);
    g = m + (A'*p) + (A'*P*c);
    h = (0.5*trace(Sig*C'*P*C)) + (0.5*c'*P*c) + (p'*c) + d + r;
    // Recursion de Riccatti
    Pnew = G - (E'*invF*E);
    pnew = g - (E'*invF*f);
    rnew = h - (0.5*f'*invF*f);
    // Controle optimal
    Unew = - (invF*E);
    unew = - (invF*f);
    //
    // La fonction de Bellman est de la forme
    //     0.5 x Pnew x  +  pnew x  +  rnew
    // et le feedback optimal est de la forme
    //     Unew x  +  unew

endfunction
