%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   		Fixed-step fourth-order Runge-Kutta method%                              
%
%                          January, 8th 2023
%       Author: Maribel Tello-Bello / Pedro Pancóatl-Bortolotti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,x,y]=RK4(ODEl,ODE2,x0,y0,h,a,b,Signal,W)
t(1) = a; 
x(1) = x0;  
y(1) = y0;
%n = (b - a)/h;
  n = length(Signal);
      for i = 1:n 
          t(i+1) = t(i)  + h; 
          tm = t(i)  + h/2; 
          Kx1 = ODEl(t(i) ,x(i) ,y(i)); 
          Ky1 = ODE2(t(i) ,x(i),y(i)) +                       (W)*Signal(i)  ; 
          
          Kx2 = ODEl(tm,x(i)+ Kx1*h/2,y(i)+  Ky1*h/2); 
          Ky2 = ODE2(tm,x(i)+ Kx1*h/2,y(i)+  Ky1*h/2) +       (W)*(Signal(i) + Ky1*Signal(i)*h/2 ); 
          
          Kx3 = ODEl(tm,x(i)+  Kx2*h/2,y(i)+  Ky2*h/2); 
          Ky3 = ODE2(tm,x(i)+ Kx2*h/2,y(i)+  Ky2*h/2) +       (W)*(Signal(i) +  Ky2*Signal(i)*h/2);
          
          Kx4 = ODEl(t(i +  1),x(i)+  Kx3*h,y(i)+  Ky3*h); 
          Ky4 = ODE2(t(i +  1),x(i)+ Kx3*h,y(i) +  Ky3*h) +   (W)*(Signal(i) + Ky3*Signal(i)*h) ;  
          
          x(i+1) = x(i)  + (Kx1 +  2*Kx2  + 2*Kx3  + Kx4)*h/6;
          if (abs(x(i + 1)) >= pi/2)
              if (x(i+1) > 0) 
                  x(i+1) = 1e-3;
              elseif x(i+1) < 0 
                  x(i+1) = -1e-3;
              end
          end

          y(i+1) = y(i)  + (Ky1 +  2*Ky2 +  2*Ky3 +  Ky4)*h/6;
          if (abs(y(i + 1)) >= pi/2)
              if (y(i+1) > 0) 
                  y(i+1) = 1e-3;
              elseif y(i+1) < 0 
                  y(i+1) = -1e-3;
              end
          end
      end 
  end
