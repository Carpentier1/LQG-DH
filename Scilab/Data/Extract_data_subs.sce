//  ----------------------------------------------------------------------
// Extraction des donnees des sous-problemes a partir des donnees globales
//  ----------------------------------------------------------------------

// On calcule a partir des donnees du probleme global celles
// de chaque sous-probleme le constituant, en extrayant les
// sous matrices propres aux sous-problemes

// Boucle sur les sous-problemes
// -----------------------------

Data_subs = list();
for ns=1:NLQG
    // Position du sous-probleme dans le probleme global
    nxf = 0;
    nuf = 0;
    for nss=1:ns
        nxf = nxf + nxss(nss);
        nuf = nuf + nuss(nss);
    end
    nxi = nxf - nxss(ns) + 1;
    nui = nuf - nuss(ns) + 1;
    // Donnees propres au sous-probleme
    //
    // - taille du sous-probleme
    nxs = nxss(ns);
    nus = nuss(ns);
    // - etat initial
    xinis = xini(nxi:nxf);
    // - cout final
    MTs = MT(nxi:nxf,nxi:nxf);
    mTs = mT(nxi:nxf);
    dTs = dT / NLQG;
    // - donnees propres aux classes de periodicite
    classess = list();
    for icl=1:Ncl
        // . recuperation des donnees de la classe
        Data = classes(icl);
        [Mt,Nt,St,mt,nt,dt,At,Bt,Ct,ct,Moy,Sig,Theta,theta] = ...
            Get_data_clas(Data);
        // . cout instantane
        Mts = Mt(nxi:nxf,nxi:nxf);
        Nts = Nt(nui:nuf,nui:nuf);
        Sts = St(nxi:nxf,nui:nuf);
        mts = mt(nxi:nxf,:);
        nts = nt(nui:nuf,:);
        dts = dt / NLQG;
        // . dynamique
        Ats = At(nxi:nxf,nxi:nxf);
        Bts = Bt(nxi:nxf,nui:nuf);
        Cts = Ct(nxi:nxf,:);
        cts = ct(nxi:nxf,:);
        // . bruit
        //Moys = Moy(nxi:nxf);
        //Sigs = Sig(nxi:nxf,nxi:nxf);
        // . terme de couplage
        Thetas = Theta(1:nc,nui:nuf);
        // . stockage des donnees de la classe
        classess($+1) = list(Mts,Nts,Sts,mts,nts,dts, ...
                             Ats,Bts,Cts,cts,Moy,Sig,Thetas,theta);
    end
    // - donnees du sous-probleme
    Data = list(nxs,nus,nc,xinis,MTs,mTs,dTs,classess);
    // - stockage des donnees de l'ensemble des sous-problemes
    Data_subs($+1) = Data;
end
