%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Detection of intermittent periods by measuring the maximums
%                    between chaos and periodicity                            
%
%                          January, 8th 2023
%       Author: Maribel Tello-Bello / Pedro Pancóatl-Bortolotti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ deltaT ] = Delta4( del )
global ts;
    
    if length(del)==1 || isempty(del)
        deltaT=[0];
    else

        for j=2:1:length(del)            %Ciclo que calcula las deltas del vector total
            deltaT(j-1)=(del(j)-del(j-1)).*ts;
        end
    
    end

end

