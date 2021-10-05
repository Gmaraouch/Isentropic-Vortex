function [Ru]=R_Central_4(F,G,dx,dy)
   
    %dw/dt=-dF/dx-dG/dy
    %dw/dt=-(-F(i+2)+8F(i+1)-8F(i-1)+F(i-2))/(2*dx)-(G(y+1)-G(y-1))/(2*dy)
    [n,m]=size(F{1}); 
    %n=number of rows (y), m=number of colums (x)
    Ru=cell(length(F),1);
    Rux=Ru;
    Ruy=Ru;
    csntx=-1/(12*dx);
    csnty=-1/(12*dy);
    
    for i=1:4
        
        for x=1:m
           
            if(x==1)
                Rux{i}(:,x)=csntx*(-F{i}(:,x+2)+8*F{i}(:,x+1)-8*F{i}(:,m)+F{i}(:,m-1));
            elseif(x==2)
                Rux{i}(:,x)=csntx*(-F{i}(:,x+2)+8*F{i}(:,x+1)-8*F{i}(:,1)+F{i}(:,m));
            elseif(x==(m-1))
                Rux{i}(:,x)=csntx*(-F{i}(:,1)+8*F{i}(:,x+1)-8*F{i}(:,x-1)+F{i}(:,x-2));
            elseif(x==m)
                Rux{i}(:,x)=csntx*(-F{i}(:,2)+8*F{i}(:,1)-8*F{i}(:,x-1)+F{i}(:,x-2));
            else
                Rux{i}(:,x)=csntx*(-F{i}(:,x+2)+8*F{i}(:,x+1)-8*F{i}(:,x-1)+F{i}(:,x-2));
            end
        end
        
        for y=1:n
           
            if(y==1)
               Ruy{i}(y,:)=csnty*(-G{i}(y+2,:)+8*G{i}(y+1,:)-8*G{i}(n,:)+G{i}(n-1,:));
            elseif(y==2)
               Ruy{i}(y,:)=csnty*(-G{i}(y+2,:)+8*G{i}(y+1,:)-8*G{i}(y-1,:)+G{i}(n,:));
            elseif(y==(n-1))
               Ruy{i}(y,:)=csnty*(-G{i}(1,:)+8*G{i}(y+1,:)-8*G{i}(y-1,:)+G{i}(y-2,:));
            elseif(y==n)
               Ruy{i}(y,:)=csnty*(-G{i}(2,:)+8*G{i}(1,:)-8*G{i}(y-1,:)+G{i}(y-2,:));
            else
               Ruy{i}(y,:)=csntx*(-G{i}(y+2,:)+8*G{i}(y+1,:)-8*G{i}(y-1,:)+G{i}(y-2,:));
            end
            
        end
        
        Ru{i}=Rux{i}+Ruy{i};
    end
end