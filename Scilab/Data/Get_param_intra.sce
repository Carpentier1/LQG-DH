function [P,Alpha,betaa,Gamma,delta,epsilon, ...
            Umat,ulin,ucst]=Get_param_intra(Data);

    // --------------------------------------------------
    // Acquisition des parametres d'une fonction intraday
    // --------------------------------------------------

    P = Data(1);
    Alpha = Data(2);
    betaa = Data(3);
    Gamma = Data(4);
    delta = Data(5);
    epsilon = Data(6);
    Umat = Data(7);
    ulin = Data(8);
    ucst = Data(9);

endfunction
