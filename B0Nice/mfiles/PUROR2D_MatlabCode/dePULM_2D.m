% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_2D mag_tmp mean_tmp] = dePULM_2D(phase_original,filter_mag, filter_th)
%----------------------------------------------------------------------
%   load the complex imaging file
    pi_2 = pi*2;
    phase_range = 4.0*pi; %[-phase_range phase_range]
    [yr xr] = size(phase_original);
    base_line = yr/2;
%----------------------------------------------------------------------
    mag_tmp = zeros(1,(yr*xr));
    for ii = 1:(yr*xr)
        %if (mag_original(ii) >= thresh_mag) && (filter_mag(ii) >= filter_th)
        if (filter_mag(ii) >= filter_th)
        mag_tmp(ii) = 1; 
        end       
    end
    %
    mag_tmp = reshape(mag_tmp,yr,xr);
    %
    %----------------------------------------------------
%    figure(4)
%    imshow(unwrapped_phase,[-pi pi])
    %pause
%---------------------------------------------------------------------
%    seg_tim = cputime - t_seg
%---------------------------------------------------------------------
    cell_signal = cell(1,yr);
    for jj = 2:(yr-1)
        sum_tmp_y = mag_tmp(jj+1,:) + mag_tmp(jj,:)+ mag_tmp(jj-1,:);
        mag_signal_tmp = zeros(1,xr);        
        kk = 1;
        for ii = 2:(xr - 1)
            if sum_tmp_y(ii) >= 2 && mag_tmp(jj,ii-1) + mag_tmp(jj,ii)+ mag_tmp(jj,ii+1) >= 2
            mag_signal_tmp(kk) = ii; 
            kk = kk + 1;
            end       
        end
        %
        if kk <= xr/4
            for ii = 2:(xr - 1)
                if sum_tmp_y(ii) >= 2
                mag_signal_tmp(kk) = ii; 
                kk = kk + 1;
                end       
            end            
        end
        mag_signal_tmp(kk:xr) = [];  
        cell_signal(1,jj) = {mag_signal_tmp};
        %
    end    
    %----------------------------------------------------------------
    phase_tmp = zeros(1,(yr*xr));
    phase_tmp = reshape(phase_tmp,yr,xr);
    for jj = 2:(yr-1)
        mag_tmp(jj,:) = 0;
        index_jj = cell_signal{1,jj};
        mag_tmp(jj,index_jj) = 1;
        phase_tmp(jj,index_jj) = phase_original(jj,index_jj);
    end
    mag_tmp(1,:) = 0; mag_tmp(yr,:) = 0; phase_tmp(1,:) = 0; phase_tmp(yr,:) = 0;
    %----------------------------------------------------------------------
    cell_connect = cell(1,yr);
    line_quality = zeros(1,yr);
    for jj = 2:(yr-1)
        sum_tmp_y = mag_tmp(jj+1,:) + mag_tmp(jj,:)+ mag_tmp(jj-1,:);
        mag_signal_tmp_connect = zeros(1,xr);
        bb = 1;
        for ii = 2:(xr - 1)
            if sum_tmp_y(ii) == 3 && mag_tmp(jj,ii-1) + mag_tmp(jj,ii)+ mag_tmp(jj,ii+1) == 3
               %sum_tmp_x = mag_tmp(jj,ii-1) + mag_tmp(jj,ii)+ mag_tmp(jj,ii+1);
               %if sum_tmp_x >= 1
                  mag_signal_tmp_connect(bb) = ii; 
                  bb = bb + 1;
               %end    
            end    
        end
        %      
        if bb <= xr/4
              for ii = 2:(xr - 1)
                if sum_tmp_y(ii) == 3
               %sum_tmp_x = mag_tmp(jj,ii-1) + mag_tmp(jj,ii)+ mag_tmp(jj,ii+1);
               %if sum_tmp_x >= 1
                  mag_signal_tmp_connect(bb) = ii; 
                  bb = bb + 1;
               %end    
                end              
              end
        end
        mag_signal_tmp_connect(bb:xr) = []; 
        %
        if bb == 1
            line_quality(jj) = -10;
        else
        line_quality(jj) = length(mag_signal_tmp_connect);
        cell_connect(1,jj) = {mag_signal_tmp_connect};
        end
    end
 %    figure(4)
 %    imshow(brg_tmp)
 %    stop
    index_signal_tmp = find(line_quality < 3);
    [Y I] = sort(diff(index_signal_tmp),'descend');
    index_y_max = index_signal_tmp(I(1)+1);
    index_y_min = index_signal_tmp(I(1));
    %index_y_min = index_signal_tmp(1) + 1;
    %index_y_max = index_signal_tmp(length(index_signal_tmp)) - 1;
    %-----------------------------------------------------------------
    if index_y_max == yr
       index_y_max = yr-1;
    end
    if index_y_min == 1
        index_y_min = 2;
    end
%--------------------------------------------------------------------
    cell_ls = cell(1,yr);
    unwrapped_phase = phase_tmp;
    mean_connect = zeros(1,yr);
    for index_y = index_y_min:index_y_max %1:length(index_signal)
        phase_x = phase_tmp(index_y,:);
        %
        index_u = cell_signal{1,index_y};
        %index_u = 1:xr;
        phase_u = phase_x(index_u);
        %phase_tmp_2 = unwrap(phase_u);
        %index_y
        [phase_tmp_1 line_ls] = dePULM_1D(phase_u,index_u);
        %index_y
        %line_ls
        %plot(phase_tmp_1,'o-r');
        %pause
        cell_ls(1,index_y) = {line_ls};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        phase_x(index_u) = phase_tmp_1(1:length(index_u)); 
        unwrapped_phase(index_y,:) = phase_x;
        %   calculating the connect_mean
        index_x = cell_connect{1,index_y};
        mean_connect(index_y) = mean_liu(phase_x(index_x));
        unwrapped_phase(index_y,index_u) = unwrapped_phase(index_y,index_u) - round(mean_connect(index_y)/pi_2)*pi_2;
        mean_connect(index_y) = mean_liu(unwrapped_phase(index_y,index_x));            
    end
    %--------------------------------------------------------
    mean_tmp = mean_connect;
    %--------------------------------------------------------
%    figure(40)
%    imshow(unwrapped_phase,[-phase_range phase_range]);
%    pause 
    %    stop
%     figure(41)
%     plot(mean_tmp/pi,'o-r'); hold off;
   % plot(mean_connect/pi,'o-g'); hold off;
   % pause
%------------------------------------------------------------------------
%    line_time = cputime - t_seg
% line unwrapping done
%--------------------------------------------------------------------------
%--------------------------------------------------------------------
%   Start line shift by unwrapping the means    
%--------------------------------------------------------------------
%   unwrap the mean values before global shfit the phase data
    mean_u = mean_tmp(index_y_min:index_y_max);
    mean_unwrap = zeros(1,yr);
    mean_unwrap(index_y_min:index_y_max) = unwrap(mean_u);
%----------------------------------------------------------
    mean_unwrap(index_y_min:index_y_max) = mean_unwrap(index_y_min:index_y_max)...
                                         - round(mean_unwrap(base_line)/pi_2)*pi_2;
%--------------------------------------------------------------------------
%   shift the phase data
    for index_y = index_y_min:index_y_max
        diff_test = mean_unwrap(index_y) - mean_tmp(index_y);
        if abs(diff_test) > pi
        unwrapped_phase(index_y,:) = unwrapped_phase(index_y,:)...
                               + pi_2*round(diff_test/(pi_2)).*mag_tmp(index_y,:);                                      
        end
        index_u = cell_connect{1,index_y};
        mean_tmp(index_y) = mean_liu(unwrapped_phase(index_y,index_u));
    end
%
%--------------------------------------------------------------------------
%   done global shift
%--------------------------------------------------------------------------
%
%   figure(41)
%    plot(mean_tmp/pi,'o-b');hold on;
    %
    figure(5)
    imshow(unwrapped_phase,[-phase_range phase_range]);
    pause
%    
%   check the quality of the bridge pixels
   quality_up = zeros(1,yr);
   quality_dw = zeros(1,yr);
    for index_y = index_y_min:index_y_max
        index_x = cell_connect{1,index_y};
        diff_up = unwrapped_phase(index_y,index_x) -  unwrapped_phase(index_y - 1,index_x);
        diff_dw = unwrapped_phase(index_y,index_x) -  unwrapped_phase(index_y + 1,index_x);
        quality_up(index_y) = length(find(abs(diff_up) < 1.0*pi))/length(index_x);
        quality_dw(index_y) = length(find(abs(diff_dw) < 1.0*pi))/length(index_x);            
    end
    quality_up(index_y_min) = 1;
    quality_dw(index_y_max) = 1;
    quality_bridge = quality_up + quality_dw;
%   
%    figure(8)
%    plot(quality_up,'o-b');hold on;
%    plot(quality_dw,'o-r');hold on;
%    plot(quality_bridge,'o-g');hold on;
%------------------------------------------------------------------------
%   find the protect region
    diff_fst_tmp = find(quality_bridge < 1.8);
    diff_fst = [index_y_min diff_fst_tmp index_y_max];
    [y protect_tmp] = sort(diff(diff_fst),'descend');
    protect_start_tmp = diff_fst(protect_tmp(1));
    protect_end_tmp = diff_fst(protect_tmp(1)+1) - 1;
    diff_2nd = find(quality_bridge(protect_start_tmp:protect_end_tmp) < 1.9) + protect_start_tmp;
    if length(diff_2nd) <= 2
    protect_start = protect_start_tmp;
    protect_end = protect_end_tmp;
    else
    [y protect_tmp] = sort(diff(diff_2nd),'descend'); 
    protect_start = diff_2nd(protect_tmp(1));
    protect_end = diff_2nd(protect_tmp(1)+1) - 1;
    end
%################################################################# 
%    index_y going down jj_loop = 1
%    index_y going up jj_loop = 2
%#################################################################
    jump_w = 1.0*pi;
    gg_w = 3;
%--------------------------------------------------------------
   for ff_loop = 1:4 
        if mod(ff_loop,2) == 1
        nxt_drn = -1;
        ref_drn = +1;
        jj = protect_start;
        else
        nxt_drn = +1;
        ref_drn = -1;
        jj = protect_end;            
        end
        do_2D = 1;
        while do_2D == 1
        %for jj = protect_start:-1:index_y_min
            if nxt_drn == -1
            trl_tmp = quality_dw(jj);
            else
            trl_tmp = quality_up(jj);    
            end
            phase_RF = unwrapped_phase(jj + ref_drn,:);
            phase_TR = unwrapped_phase(jj,:);            
            phase_NT = unwrapped_phase(jj + nxt_drn,:);
            flag_display = 0;
            %jj
            index_jj = cell_connect{1,jj};
            %-------------------------------------------------------------
            %figure(9)
            %plot(phase_RF,'o-k'); hold on;
            %plot(phase_TR,'o-r'); hold on;
            %plot(phase_NT,'o-g'); hold off;
            %pause
            %       Method 2
            tmp_tmp = zeros(1,xr);
            tmp_tmp(index_jj) = 1;
            mag_jj = mag_tmp(jj,:);
            diff_test = (phase_RF - phase_TR).*tmp_tmp;
            jump_tmp  = find(abs(diff_test) > jump_w);
            if isempty(jump_tmp) || trl_tmp > 0.95
            else
            %--------------------------------------------------------------
            [phase_jj_out] = dePU_line_merging(phase_TR,mag_jj,index_jj,diff_test,jump_tmp,gg_w,xr);
            tmp_line_phase = phase_jj_out; 
            unwrapped_phase(jj,:) = tmp_line_phase;
            end
            jj = jj + nxt_drn;
            %
            %
            if jj < protect_start_tmp
            do_2D = 0;
            end
            if jj > protect_end_tmp
            do_2D = 0;
            end
            %  
        end     
   end
   %
   for index_y = protect_start_tmp:protect_end_tmp
        index_x = cell_connect{1,index_y};
        diff_up = unwrapped_phase(index_y,index_x) -  unwrapped_phase(index_y - 1,index_x);
        diff_dw = unwrapped_phase(index_y,index_x) -  unwrapped_phase(index_y + 1,index_x);
        quality_up(index_y) = length(find(abs(diff_up) < 1.0*pi))/length(index_x);
        quality_dw(index_y) = length(find(abs(diff_dw) < 1.0*pi))/length(index_x);
        quality_bridge(index_y) = quality_up(index_y) + quality_dw(index_y);
    end
%--------------------------------------------------------------
   for jj_loop = 1:2 
        if mod(jj_loop,2) == 1
        nxt_drn = -1;
        ref_drn = +1;
        jj = protect_start_tmp;
        else
        nxt_drn = +1;
        ref_drn = -1;
        jj = protect_end_tmp;            
        end
        do_2D = 1;
        while do_2D == 1
        %for jj = protect_start:-1:index_y_min
            phase_RF = unwrapped_phase(jj + ref_drn,:);
            phase_TR = unwrapped_phase(jj,:);            
            flag_display = 0;
            %jj
            index_jj = cell_connect{1,jj};
            ref_phase = unwrapped_phase(jj+ref_drn,index_jj);
            trl_phase = unwrapped_phase(jj,index_jj);
            nxt_phase = unwrapped_phase(jj + nxt_drn,index_jj);
            diff_nxt = trl_phase - nxt_phase;
            diff_ref = trl_phase - ref_phase;
            diff_nxt_ref = nxt_phase - ref_phase;
            trl_nxt = length(find(abs(diff_nxt) < 1.0*pi))/length(index_jj);
            trl_ref = length(find(abs(diff_ref) < 1.0*pi))/length(index_jj);
            nxt_ref = length(find(abs(diff_nxt_ref) < 1.0*pi))/length(index_jj);            
            if mod(jj_loop,2) == 1
            quality_up(jj) = trl_nxt;
            quality_dw(jj) = trl_ref;
            else
            quality_up(jj) = trl_ref;
            quality_dw(jj) = trl_nxt;                
            end
            quality_bridge(jj) = trl_nxt + trl_ref;
            if jj_loop <= 2
                th_trl = 0.9;
            else
                th_trl = 0.8;
            end
            %
        if  trl_ref <= th_trl 
            if flag_display == 1
            display ('need 2D');
            display ('step 0');
                trl_nxt
                trl_ref
                nxt_ref
            figure(9)
            plot(ref_phase,'o-k'); hold on;
            plot(trl_phase,'o-r'); hold on;
            plot(nxt_phase,'o-g'); hold off;
            pause
            end
          if jj_loop <= 2 
            %--------------------------------------------------------------
            ls_mag = cell_signal{1,jj};
            %--------------------------------------------------------------
            ls_phase = phase_TR(ls_mag);
            re_phase = phase_RF(ls_mag);
            [out_TR] = dePULM_line_merging(ls_phase,cell_ls{1,jj},re_phase);
            tmp_line_phase = zeros(1,xr);
            tmp_line_phase(ls_mag) = out_TR;
            %--------------------------------------------------------------
            %recaculated the quality
            trl_phase = tmp_line_phase(index_jj);
            diff_nxt = trl_phase - nxt_phase;
            diff_ref = trl_phase - ref_phase;
            diff_nxt_ref = nxt_phase - ref_phase;
            trl_nxt = length(find(abs(diff_nxt) < 1.0*pi))/length(index_jj);
            trl_ref = length(find(abs(diff_ref) < 1.0*pi))/length(index_jj);
            nxt_ref = length(find(abs(diff_nxt_ref) < 1.0*pi))/length(index_jj);
            if mod(jj_loop,2) == 1
            quality_up(jj) = trl_nxt;
            quality_dw(jj) = trl_ref;
            else
            quality_up(jj) = trl_ref;
            quality_dw(jj) = trl_nxt;                
            end
            quality_bridge(jj) = trl_nxt + trl_ref;
            %--------------------------------------------------------------
            method_1_out = tmp_line_phase;
            method_1_trl_ref = trl_ref;
            %--------------------------------------------------------------
                if flag_display == 1
                display ('method 1');
                trl_nxt
                trl_ref
                nxt_ref
                figure(10)
                plot(ref_phase,'o-k'); hold on;
                plot(trl_phase,'o-r'); hold on;
                plot(nxt_phase,'o-g'); hold off;           
                pause
                end
            %-------------------------------------------------------------
          else
            %       Method 2
            tmp_tmp = zeros(1,xr);
            tmp_tmp(index_jj) = 1;
            mag_jj = mag_tmp(jj,:);
            diff_test = (phase_RF - phase_TR).*tmp_tmp;
            jump_tmp  = find(abs(diff_test) > jump_w);
            %--------------------------------------------------------------
            if isempty(jump_tmp)
            tmp_line_phase = phase_TR;
            else
            [phase_jj_out] = dePU_line_merging(phase_TR,mag_jj,index_jj,diff_test,jump_tmp,gg_w,xr);
            tmp_line_phase = phase_jj_out;
            end
            %
            trl_phase = tmp_line_phase(index_jj);
            diff_nxt = trl_phase - nxt_phase;
            diff_ref = trl_phase - ref_phase;
            diff_nxt_ref = nxt_phase - ref_phase;
            trl_nxt = length(find(abs(diff_nxt) < 1.0*pi))/length(index_jj);
            trl_ref = length(find(abs(diff_ref) < 1.0*pi))/length(index_jj);
            nxt_ref = length(find(abs(diff_nxt_ref) < 1.0*pi))/length(index_jj);            
            if mod(jj_loop,2) == 1
            quality_up(jj) = trl_nxt;
            quality_dw(jj) = trl_ref;
            else
            quality_up(jj) = trl_ref;
            quality_dw(jj) = trl_nxt;                
            end
            quality_bridge(jj) = trl_nxt + trl_ref;
                if flag_display == 1
                display ('method 2');
                trl_nxt
                trl_ref
                nxt_ref
                figure(11)
                plot(ref_phase,'o-k'); hold on;
                plot(trl_phase,'o-r'); hold on;
                plot(nxt_phase,'o-g'); hold off;
                pause
                end
            %--------------------------------------------------------------
            method_2_out = tmp_line_phase;
          end
            unwrapped_phase(jj,:) = tmp_line_phase;
        end
            %
            %index_jj = cell_connect{1,jj};
            jj_tmp = unwrapped_phase(jj,index_jj);
            mean_tmp(jj) = mean_liu(jj_tmp);  
            %
            if mod(jj_loop,2) ==1
                if jj > (index_y_min + 2)
                mean_diff = mean_tmp(jj + ref_drn) - mean_tmp(jj);
                    if abs(mean_diff) > jump_w
                    kk = jj;
                    unwrapped_phase(index_y_min:kk,:) = unwrapped_phase(index_y_min:kk,:)...
                               + pi_2*round(mean_diff/(pi_2)).*mag_tmp(index_y_min:kk,:);
                    mean_tmp(index_y_min:kk) = mean_tmp(index_y_min:kk)...
                               + pi_2*round(mean_diff/(pi_2));
                    end
                end
                            %
               if jj > (index_y_min + 2)
               mean_diff = mean_tmp(jj) - mean_tmp(jj +nxt_drn);
                    if abs(mean_diff) > jump_w
                    kk = jj + nxt_drn;
                    unwrapped_phase(index_y_min:kk,:) = unwrapped_phase(index_y_min:kk,:)...
                               + pi_2*round(mean_diff/(pi_2)).*mag_tmp(index_y_min:kk,:);
                    mean_tmp(index_y_min:kk) = mean_tmp(index_y_min:kk)...
                               + pi_2*round(mean_diff/(pi_2));
                    end
               end
               %
            else
                %
                if jj < (index_y_max - 2)
                mean_diff = mean_tmp(jj + ref_drn) - mean_tmp(jj);
                    if abs(mean_diff) > jump_w*1.2
                    kk = jj;
                unwrapped_phase(kk:index_y_max,:) = unwrapped_phase(kk:index_y_max,:)...
                               + pi_2*round(mean_diff/(pi_2)).*mag_tmp(kk:index_y_max,:);
                mean_tmp(kk:index_y_max) = mean_tmp(kk:index_y_max)...
                               + pi_2*round(mean_diff/(pi_2));
                    end
                end
                %
                if jj < (index_y_max - 2)
                mean_diff = mean_tmp(jj) - mean_tmp(jj + nxt_drn);
                    if abs(mean_diff) > jump_w*1.2
                    kk = jj + nxt_drn;
                    unwrapped_phase(kk:index_y_max,:) = unwrapped_phase(kk:index_y_max,:)...
                               + pi_2*round(mean_diff/(pi_2)).*mag_tmp(kk:index_y_max,:);
                    mean_tmp(kk:index_y_max) = mean_tmp(kk:index_y_max)...
                               + pi_2*round(mean_diff/(pi_2));
                    end
                end
            end
            %--------------------------------------------------------------
            jj = jj + nxt_drn;
            if jj < index_y_min
            do_2D = 0;
            end
            if jj > index_y_max
            do_2D = 0;
            end
            %  
        end     
   end       
%    figure(41)
%    plot(mean_tmp/pi,'o-r');  hold off
%    plot(diff(mean_tmp)/pi,'o-k');  hold off;
    %-------------------------------------------------
%--------------------------------------------------------------------
%   For the final step, Start line shift by unwrapping the means    
%--------------------------------------------------------------------
    quality_bridge = zeros(1,yr);
    for index_y = (index_y_min+1):(index_y_max-1)
        index_x = cell_connect{1,index_y};
        diff_up = unwrapped_phase(index_y,index_x) -  unwrapped_phase(index_y - 1,index_x);
        diff_dw = unwrapped_phase(index_y,index_x) -  unwrapped_phase(index_y + 1,index_x);
        quality_bridge(index_y) = (length(find(abs(diff_up) < 1.0*pi)) + length(find(abs(diff_dw) < 1.0*pi)))/length(index_x);            
    end
%   
    figure(8)
    plot(quality_bridge,'o-b');hold off;
    %
out_2D = unwrapped_phase;
%--------------------------------------------------------------------------
%   phase_unwrapping was done
