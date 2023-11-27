function [newtime] = f_newtime_NREL(time_vector_UTC,round_sec,reverse)
%function [newtime] = f_newtime_NREL(time_vector_UTC,round_sec,reverse)
% Format Time Vector into columns to meet format specs by NREL
% 
% Inputs
%       time_vector_UTC =  A datetime array
%             round_sec =  1 if yes round to nearest second | 0 for no
%                          rounding
%             reverse   =  1 of you want to convert back to one cell timestamp 
%                          0 if you want to go to 6 cell timestamp (NREL
%                          Format)
%
% Output
%       newtime =  A datetime array
%
%       YYYY | MM | dd | HH | mm | ss.sssssssss 

if reverse == 0

    % Seperate out years
    newtime(:,1) = year(time_vector_UTC);
    % Seperate out months
    newtime(:,2) = month(time_vector_UTC);
    % Seperate out Days
    newtime(:,3) = day(time_vector_UTC);
    % Seperate out Hours
    newtime(:,4) = hour(time_vector_UTC);
    % Seperate out Minutes
    newtime(:,5) = minute(time_vector_UTC);   
    % Seperate out Seconds
    if round_sec == 1
        newtime(:,6) = round(second(time_vector_UTC));

    elseif round_sec == 0
        newtime(:,6) = second(time_vector_UTC);
        
    else
        return
    end
    
elseif reverse == 1
    
    newtime = datetime(time_vector_UTC(:,1),time_vector_UTC(:,2),...
        time_vector_UTC(:,3),time_vector_UTC(:,4),time_vector_UTC(:,5),...
        time_vector_UTC(:,6),'format','MM/dd/yyyy HH:mm:ss.SSSSSSSSS');
end
end

