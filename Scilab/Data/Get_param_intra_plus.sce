function [E,F,G,Alphah,betaah]=Get_param_intra_plus(Data);

    // -----------------------------------------------------------------
    // Acquisition de parametres supplementaires d'une fonction intraday
    // dans le cas ou la recursion de Riccatti ne tient pas compte de la
    // contrainte couplante
    // -----------------------------------------------------------------

    E = Data(1);
    F = Data(2);
    G = Data(3);
    Alphah = Data(4);
    betaah = Data(5);

endfunction
