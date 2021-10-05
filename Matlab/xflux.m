function F=xflux(w,Ma,gamma)

%w=[rho, rho*u,rho*v,rho*E];
%F=[rho*u,rho*u^2+p,rho*u*v,u*(rho*E+p)];

rho=w{1};
rhou=w{2};
rhov=w{3};
rhoE=w{4};

u=rhou./rho;
v=rhov./rho;

%p=(rho.^(gamma))/(gamma*Ma^2);
p=(gamma-1)*(rhoE-(1/2).*rho.*(u.^2+v.^2));
F=cell(4,1);

F{1}=rhou;
F{2}=rhou.*u+p;
F{3}=rhou.*v;
F{4}=u.*(rhoE+p);

end