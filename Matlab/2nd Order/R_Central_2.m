function [Ru]=R_Central_2(F,G,dx,dy)
   
    %dw/dt=-dF/dx-dG/dy
    %dw/dt=-(F(i+1)-F(i-1)/(2*dx)-(G(y+1)-G(y-1))/(2*dy)
    [n,m]=size(F{1}); 
    %n=number of rows (y), m=number of colums (x)
    Ru=cell(length(F),1);
    Rux=Ru;
    Ruy=Ru;
    csntx=-1/(2*dx);
    csnty=-1/(2*dy);
    
    for i=1:4
        
        for x=1:m
           
            if(x==1)
                Rux{i}(:,x)=csntx*(F{i}(:,x+1)-F{i}(:,m));
            elseif(x==m)
                Rux{i}(:,x)=csntx*(F{i}(:,1)-F{i}(:,x-1));
            else
                Rux{i}(:,x)=csntx*(F{i}(:,x+1)-F{i}(:,x-1));
            end
        end
        
        for y=1:n
           
            if(y==1)
               Ruy{i}(y,:)=csnty*(G{i}(y+1,:)-G{i}(n,:));
            elseif(y==n)
               Ruy{i}(y,:)=csnty*(G{i}(1,:)-G{i}(y-1,:));
            else
               Ruy{i}(y,:)=csnty*(G{i}(y+1,:)-G{i}(y-1,:));
            end
            
        end
        
        Ru{i}=Rux{i}+Ruy{i};
    end
end