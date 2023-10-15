%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Detection of intermittent periods by measuring the change of transitions
%                       between chaos and periodicity                            
%
%                          January, 8th 2023
%       Author: Maribel Tello-Bello / Pedro Pancóatl-Bortolotti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ deltaT ] = Delta( umbral, vec )

global ts
n=0;                    
tam=size(vec);          
for i=1:1:tam           
    if vec(i)>umbral
        n=n+1;          
        tiempos(n)=i;   
    end;      
end;
x=0;m=1;k=1;            
tamti=length(tiempos);  

for i=2:1:tamti         
    dif=tiempos(i)-tiempos(i-1);
    if dif<600          
       x=x+1;           
    else                
       m=round(x/2);    
       delta(k)= tiempos(i-m);  
       x=0; k=k+1;     
    end;
 end;
m=round(x/2); delta(k)= tiempos(i-m);   
if k>=2
    for j=2:1:k             
        deltaT(j-1)=(delta(j)-delta(j-1)).*ts;
    end
else
   
   deltaT=[0];
end
end

