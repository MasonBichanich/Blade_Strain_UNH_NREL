function [Total_pts_flagged,MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct)

%Master Flag Matrix
    MasterFlag = ones(length(Flagged_Struct.SIGSTR),width(Flagged_Struct.SIGSTR));

% for i = 1:length(Flagged_Struct.SIGSTR)
%     
% 
%         if Flagged_Struct.SIGSTR(i) == 0 || Flagged_Struct.Corr(i) == 0 || Flagged_Struct.Horz_Current_Mag(i) == 0 ||...
%            Flagged_Struct.Horz_Current_Dir(i) == 0 || Flagged_Struct.INSVELX_max(i) == 0 ||  Flagged_Struct.INSVELY_max(i) == 0 ||...
%             Flagged_Struct.INSVELZ_max(i) == 0 || Flagged_Struct.INSVELX_dt(i) == 0 || Flagged_Struct.INSVELY_dt(i) == 0 || ...
%             Flagged_Struct.INSVELX_spk(i) == 0 || Flagged_Struct.INSVELY_spk(i) == 0 || Flagged_Struct.INSVEL_stdfltr(i) == 0
%              
%                        MasterFlag(i,1) = 0;
%          
%         else
%             MasterFlag(i,1) = 1;
%         end
%     
% end
 
MasterFlag( Flagged_Struct.SIGSTR == 0 | Flagged_Struct.Corr == 0 | Flagged_Struct.Horz_Current_Mag == 0 |...
       Flagged_Struct.Horz_Current_Dir == 0 | Flagged_Struct.INSVELX_max == 0 |  Flagged_Struct.INSVELY_max == 0 |...
          Flagged_Struct.INSVELZ_max == 0 | Flagged_Struct.INSVELX_dt == 0 | Flagged_Struct.INSVELY_dt == 0 | ...
            Flagged_Struct.INSVELX_spk == 0 | Flagged_Struct.INSVELY_spk == 0 | Flagged_Struct.INSVEL_stdfltr == 0) = 0;
              

Master_Flagged = find(MasterFlag(:,1)==0);
 
        
Total_pts_flagged = length(Master_Flagged);

end

