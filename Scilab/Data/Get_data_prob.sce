function [nx,nu,nc,xini,MT,mT,dT,classes]=Get_data_prob(Data);

    // -----------------------------------
    // Acquisition des donnees du probleme
    // -----------------------------------

    // Cette acquisition ne porte que sur les donnees susceptibles
    // d'etre modifiees, par exemple dans un sous-probleme lors de
    // la decomposition en espace

    nx = Data(1);
    nu = Data(2);
    nc = Data(3);
    xini = Data(4);
    MT = Data(5);
    mT = Data(6);
    dT = Data(7);
    classes = Data(8);

endfunction
