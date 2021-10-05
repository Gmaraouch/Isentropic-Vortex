function [wt,F,G]=fRK44_2(w,dx,dy,dt,Ma,gamma)
    
    w2=cell(4,1);
    w3=w2;
    w4=w2;
    wt=w2;
    
    w1=w;
    F1=xflux(w1,Ma,gamma);
    G1=yflux(w1,Ma,gamma);
    Rw1=R_Central_2(F1,G1,dx,dy);
    
    for i=1:4
    w2{i}=w{i}+(dt/2).*Rw1{i};
    end
    
    F2=xflux(w2,Ma,gamma);
    G2=yflux(w2,Ma,gamma);
    Rw2=R_Central_2(F2,G2,dx,dy);
    
    for i=1:4
    w3{i}=w{i}+(dt/2).*Rw2{i};
    end
    
    F3=xflux(w3,Ma,gamma);
    G3=yflux(w3,Ma,gamma);
    Rw3=R_Central_2(F3,G3,dx,dy);
    
    for i=1:4
    w4{i}=w{i}+(dt).*Rw3{i};
    end
    
    F4=xflux(w4,Ma,gamma);
    G4=yflux(w4,Ma,gamma);
    Rw4=R_Central_2(F4,G4,dx,dy);
    
    for i=1:4
    wt{i}=w{i}+(dt/6).*(Rw1{i}+2.*Rw2{i}+2.*Rw3{i}+Rw4{i});
    end
    
    F=xflux(w,Ma,gamma);
    G=yflux(w,Ma,gamma);
    
end