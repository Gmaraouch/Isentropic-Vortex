clear 
clc
clf

addpath('C:\Users\Gmara\OneDrive - Concordia University - Canada\Concordia\Fall 2019\ENGR 6251\Assignment 3\Question 3');

Lx=20;
Ly=20;

gamma=1.4;
R=1.5;
Ma=0.4;
xc=10;
yc=10;
S=13.5;

Nx=[60,120,240,480];
Ny=Nx;
t=20;
dt=0.01;
nt=t/dt;

Error(1:length(Nx))=0;
for i=1:length(Nx)
    dx=Lx/Nx(i);
    dy=Ly/Ny(i);
    
    x=0:dx:(Lx-dx);
    y=(0:dy:(Ly-dy));
    
    [x,y]=meshgrid(x,y);
    
    f=(1-(x-xc).^2-(y-yc).^2)/(2*R^2);
    rho=(1-((S^2*Ma^2*(gamma-1)).*exp(2*f))/(8*pi^2)).^(1/(gamma-1));
    u=S*(y-yc).*exp(f)/(2*pi*R);
    v=1-S*(x-xc).*exp(f)/(2*pi*R);
    p=(rho.^(gamma))./(gamma*Ma^2);
    
    rhoE=p./(gamma-1)+1/2*rho.*(u.^2+v.^2);
    
    w=cell(4,1);
    w{1}=rho;
    w{2}=rho.*u;
    w{3}=rho.*v;
    w{4}=rhoE;
    %{
    figure(1)
    surf(x,y,w{1});
    figure(2)
    surf(x,y,w{2});
    figure(3)
    surf(x,y,w{3});
    figure(4)
    surf(x,y,w{4});
    %}
    F=xflux(w,Ma,gamma);
    G=yflux(w,Ma,gamma);
    
    wo=w;
    %{
    wt=cell(nt,1);
    wt{1}.w=w;
    wt{1}.F=F;
    wt{1}.G=G;
 
    for t=1:nt
        
        [wt{t+1}.w,wt{t+1}.F,wt{t+1}.G]=fRK44_2(w,dx,dy,dt,Ma,gamma);
        drawnow;
        surf(x,y,wt{t+1}.w{1});
        w=wt{t+1}.w;
    end
    
    Error(i)=sqrt(sum(sum((wt{nt+1}.w{1}-wt{1}.w{1}).^2))/(Nx(i)*Ny(i)));
    %}
    
    for t=1:nt
        
        [wt]=fRK44_2(w,dx,dy,dt,Ma,gamma);
        drawnow;
        surf(x,y,wt{1});
        w=wt;
    end
    
    Error(i)=sqrt(sum(sum((w{1}-wo{1}).^2))/(Nx(i)*Ny(i)));
    
    figure(1)
    subplot(2,1,1)
    contourf(x,y,wo{1});
    title({'Isentropic Vortex Initial Solution',['With Nx=' num2str(Nx(i)) ' and Ny=' num2str(Ny(i))]});
    colorbar;
    subplot(2,1,2)
    contourf(x,y,w{1});
    colorbar;
    xlabel('x(units)');
    ylabel('y(units)');
    title({'Isentropic Vortex Using Second Order Central Spatial',['With Nx=' num2str(Nx(i)) ' and Ny=' num2str(Ny(i))]});
    saveas(gcf,['Contour_Vortex_2_Nx=' num2str(Nx(i)) '.jpg']);
    
    subplot(2,2,1)
    contourf(x,y,wo{1});
    colorbar;
    
    xlabel('x(units)');
    ylabel('y(units)');
    
    figure(2)
    surf(x,y,w{1},'EdgeColor','None');
    colorbar;
    xlabel('x(units)');
    ylabel('y(units)');
    title({'Isentropic Vortex Using Second Order Central Spatial',['With Nx=' num2str(Nx(i)) ' and Ny=' num2str(Ny(i))]});
    saveas(gcf,['Surf_Vortex_2_Nx=' num2str(Nx(i)) '.jpg']);
    clf
    
end

loglog(Nx,Error);
xlabel('Number of Grid points');
ylabel('Error');
title({'Loglog Plot of the Error Versus Number of Grid Points', 'Using Second Order Central Spatial'});
saveas(gcf,['Error_Vortex_2_Nx=' num2str(Nx(i)) '.jpg']);

slope(1:(length(Nx)-1))=(log(Error(2:length(Nx)))-log(Error(1:(length(Nx)-1))))/(log(Nx(2:length(Nx)))-log(Nx(1:(length(Nx)-1))));
