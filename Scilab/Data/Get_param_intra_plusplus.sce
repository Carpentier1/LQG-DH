function [E,F,G,K,Aphah,betaah,Alphap,betaap]=Get_param_intra_plusplus(Data)


    // -----------------------------------------------------------------
    // Acquisition de parametres supplementaires d'une fonction intraday
    // lorsque la recursion de Riccatti prend en compte la contrainte
    // couplante
    // -----------------------------------------------------------------

    E = Data(1);
    F = Data(2);
    G = Data(3);
    K = Data(4);
    Alphah = Data(5);
    betaah = Data(6);
    Alphap = Data(7);
    betaap = Data(8);

endfunction
