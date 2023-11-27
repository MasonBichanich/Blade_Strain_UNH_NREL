function [Current_Magnitude,Current_Direction,heading,theta_1r,East,North] = f_ADV_Speed_Direction_Pat(ADV_VEL)
%% The purpose of this function is to take ADV instream x and y velocities
% that have been converted to NOAA sign convention and convert them 
% into a velocity magnitude and direction

    INS_X_VEL = ADV_VEL(:,1);

    INS_Y_VEL = ADV_VEL(:,2);

%Relate x-velocity and y-velocity to east-velocity and north-velocity 

heading=105.05; %deg Angle of bridge pier, it is assumed that the platform is parallel with the pier and that the ADV is straight in the mount

theta_d = heading-90;
theta_r = deg2rad(theta_d);

North =  INS_X_VEL.*sin(theta_r) - (INS_Y_VEL.*cos(theta_r));
East  = -1.*INS_X_VEL.*cos(theta_r) - (INS_Y_VEL.*sin(theta_r));

theta_1r = atan2(North,East);

r = sqrt((East.^2) + (North.^2));

Speed_Direction(:,1)=r; %Speed

% Convert to North Heading Reference frame
Speed_Direction(:,2)=90 - rad2deg(theta_1r); %Direction

 
%Make all directions positive
for i=1:length(Speed_Direction(:,2))
        if Speed_Direction(i,2) < 0
            Speed_Direction(i,2) = 360 + Speed_Direction(i,2);
        else
        end
end

Current_Magnitude = Speed_Direction(:,1);
Current_Direction = Speed_Direction(:,2);
end


%%


% ADV Code

% quad = 1000.*ones(length(INS_X_VEL),1);
% theta = 1000.*ones(length(INS_X_VEL),1);
% rel_theta = 1000.*ones(length(INS_X_VEL),1);
% 
% R = zeros(length(INS_X_VEL),1);
% 
% for i =1:length(INS_X_VEL)
% 
%     if sign(INS_X_VEL(i)) == 0 && sign(INS_Y_VEL(i)) ~= 0
%         
%         if sign(INS_Y_VEL(i)) == 1
%             
%             rel_theta(i) = 0;
%             R(i) = abs(INS_Y_VEL(i));
%             quad(i) = 1;
%             theta(i) = heading + 90 ;    
%             
%         elseif sign(INS_Y_VEL(i)) == -1
%             
%             rel_theta(i) = 0;
%             R(i) = abs(INS_Y_VEL(i));
%             quad(i) = 2;
%             theta(i) = heading - 90 ;
%             
%         end        
%     elseif sign(INS_X_VEL(i)) ~= 0 && sign(INS_Y_VEL(i)) == 0
%         
%         if sign(INS_X_VEL(i)) == 1
%             
%             rel_theta(i) = 0;
%             R(i) = abs(INS_X_VEL(i));
%             quad(i) = 2;
%             theta(i) = heading + 180 ;    
%             
%         elseif sign(INS_X_VEL(i)) == -1
%             
%             rel_theta(i) = 0;
%             R(i) = abs(INS_X_VEL(i));
%             quad(i) = 2;
%             theta(i) = heading ;
%             
%         end
%         
%     elseif sign(INS_X_VEL(i)) == 0 && sign(INS_Y_VEL(i)) == 0
% 
%         rel_theta(i) = 0;
%         R(i) = 0;
%         quad(i) = 0;
%         theta(i) = heading ;
%         
%     else
%         
%         
%         rel_theta(i) = atand(abs(INS_Y_VEL(i)/INS_X_VEL(i)));
%         R(i) = sqrt( (INS_X_VEL(i)^2) + (INS_Y_VEL(i)^2));
% 
%         if sign(INS_X_VEL(i)) == 1 && sign(INS_Y_VEL(i)) == 1
% 
%             quad(i) = 1;
%             theta(i) = heading + 180 - rel_theta(i);
% 
%         elseif sign(INS_X_VEL(i)) == 1 && sign(INS_Y_VEL(i)) == -1
% 
%             quad(i) = 4;
%             theta(i) = heading + 180 + rel_theta(i);
% 
%         elseif sign(INS_X_VEL(i)) == -1 && sign(INS_Y_VEL(i)) == 1
% 
%             quad(i) = 2;
%             theta(i) = heading + rel_theta(i) ;
% 
%         elseif sign(INS_X_VEL(i)) == -1 && sign(INS_Y_VEL(i)) == -1
% 
%             quad(i) = 3;
%             theta(i) = heading - rel_theta(i);
% 
%         end
%     end
% end
% 
%     
% Speed_Direction(:,1)=R; %Speed
% Speed_Direction(:,2)=theta;  %Direction
    
