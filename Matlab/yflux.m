function G=yflux(w,Ma,gamma)

%w=[rho, rho*u,rho*v,rho*E];
%G=[rho*v,rho*u*v,rho*v^2+p,v*(rho*E+p)];

rho=w{1};
rhou=w{2};
rhov=w{3};
rhoE=w{4};

u=rhou./rho;
v=rhov./rho;

%p=(rho.^(gamma))/(gamma*Ma^2);
p=(gamma-1)*(rhoE-1/2*rho.*(u.^2+v.^2));
G=cell(4,1);

G{1}=rhov;
G{2}=rhov.*u;
G{3}=rhov.*v+p;
G{4}=v.*(rhoE+p);

end