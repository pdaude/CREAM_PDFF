
function [unwrap_phase] = PUROR2D_B0Matlab(data_RawB0,mask4unwrap,mask4supp,mask4STAS)
%
%----------------------------------------------------------------------
%   2D unwrapping with PUROR 
%-----------------------------------------------------------------------
[yr, xr, s] = size(data_RawB0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pi_2 = 2*pi;
%----------------------------------------------------------------------
%   Start phase unwrapping
%----------------------------------------------------------------------
slice_start = 1;
slice_end = s;
%
unwrap_phase = zeros(yr, xr, (slice_end-slice_start+1));
%
for index_slice = slice_start:slice_end   
%  
data_slice = data_RawB0(:,:,index_slice); 
%------------------------------
MASK_2D_unwrap = mask4unwrap(:,:,index_slice);
MASK_2D_supp = mask4supp(:,:,index_slice);
MASK_stas = mask4STAS(:,:,index_slice);
%-----------------------------------------------------------------
[unwrapped_phase_x, unwrapped_phase_y, mean_y, mean_x] = dePULM_2D_doit(angle(data_slice),MASK_2D_unwrap,MASK_2D_supp);
%--------------------------------------------------------------------------
phase_tmp_mean = unwrapped_phase_x(:);
phase_tmp_mean(MASK_stas(:) == 0) = [];
phase_tmp_mean(isnan(phase_tmp_mean)) = [];
tmp_mean = mean(phase_tmp_mean);
unwrapped_phase_x = unwrapped_phase_x - pi_2*round(tmp_mean/pi_2);

    unwrap_phase(:,:,index_slice) = unwrapped_phase_x;  

end % for index_slice 
%-------------------------------------
    %
end
%---------------------------------------------------------------------
