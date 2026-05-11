function [Dscis,TDscis,Nclscis,Aclscis,TCscis,TDiscis,TDfscis, ...
                classesscis,Data_globscis]= ...
                            Data_global_scission(D,TD,Ncl,Acl,TC,TDi,TDf, ...
                                                 classes,coefscis)

// ---------------------------------------------------------------
// Donnees du probleme global avec des blocs temporels plus petits
// ---------------------------------------------------------------

// Nombre de blocs apres scission
Dscis = D * coefscis;
// Correspondance entre anciens et nouveaux blocs
corbloc = zeros(D,coefscis);
dscis = 0;
for d=1:D
    for cscis=1:coefscis
        dscis = dscis + 1;
        corbloc(d,cscis) = dscis;
    end
end
// Duree des blocs apres scission
TDscis = zeros(1,Dscis);
dscis = 0;
for d=1:D
    for cscis=1:coefscis
        dscis = dscis + 1;
        TDscis(dscis) = TD(d) / coefscis;
    end
end
// Nombre de classes de periodicite apres scission
Nclscis = Ncl * coefscis;
// Correspondance entre anciennes et nouvelles classes
corclas = zeros(Ncl,coefscis);
corclasinv = zeros(Nclscis,2)
iclscis = 0;
for icl=1:Ncl
    for cscis=1:coefscis
        iclscis = iclscis + 1;
        corclas(icl,cscis) = iclscis;
        corclasinv(iclscis,1) = icl;
        corclasinv(iclscis,2) = cscis;
    end
end
// Affectation des blocs aux classes apres scission
Aclscis = zeros(1,Dscis);
dscis = 0;
for d=1:D
    icl = Acl(d);
    for cscis=1:coefscis
        dscis = dscis + 1;
        Aclscis(dscis) = corclas(icl,cscis);
    end
end
// Duree des blocs d'une meme classe apres scission
TCscis = zeros(1,Nclscis);
for iclscis=1:Nclscis
    Tclscis = 0;
    for dscis=1:Dscis
        if Aclscis(dscis)==iclscis then
            if Tclscis==0 then
                Tclscis = TDscis(dscis);
            else
                if Tclscis~=TDscis(dscis) then
                    disp('Data_Global : durees de bloc incompatibles!');
                    pause
                end
            end
        end
    end
    TCscis(iclscis) = Tclscis;
end
// Instants initial et final des blocs apres scission
TDiscis = zeros(Dscis,1);
TDfscis = zeros(Dscis,1);
t = 0;
for dscis=1:Dscis
    TDiscis(dscis) = t +1;
    TDfscis(dscis) = t + TDscis(dscis);
    t = t + TDscis(dscis);
end
// Classes de periodicite apres scission
classesscis = list();
for iclscis=1:Nclscis
    icl = corclasinv(iclscis,1);
    Data = classes(icl);
    [Mt,Nt,St,mt,nt,dt,At,Bt,Ct,ct,Moy,Sig,Theta,theta]= Get_data_clas(Data);
    cscis = corclasinv(iclscis,2);
    tiscis = ((cscis-1)*TCscis(iclscis)) + 1;
    tfscis = cscis * TCscis(iclscis);
    mtscis = mt(:,tiscis:tfscis);
    ntscis = nt(:,tiscis:tfscis);
    dtscis = dt(:,tiscis:tfscis);
    ctscis = ct(:,tiscis:tfscis);
    thetascis = theta(:,tiscis:tfscis);
    classesscis($+1) = list(Mt,Nt,St,mtscis,ntscis,dtscis, ...
                            At,Bt,Ct,ctscis,Moy,Sig,Theta,thetascis);
end

// Stockage des donnees du probleme global apres scission
Data_globscis = list(nx,nu,nc,xini,MT,mT,dT,classesscis);
endfunction
