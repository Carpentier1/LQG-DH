function [Mt,Nt,St,mt,nt,dt,At,Bt,Ct,ct,Moy,Sig,Theta,theta]= ...
             Get_data_clas(Data);

    // ---------------------------------------------------
    // Acquisition des donnees d'une classe de periodicite
    // ---------------------------------------------------

    Mt = Data(1);
    Nt = Data(2);
    St = Data(3);
    mt = Data(4);
    nt = Data(5);
    dt = Data(6);
    At = Data(7);
    Bt = Data(8);
    Ct = Data(9);
    ct = Data(10);
    Moy = Data(11);
    Sig = Data(12);
    Theta = Data(13);
    theta = Data(14);

endfunction
