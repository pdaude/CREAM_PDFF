function [line_tmp_phase] = dePULM_line_merging(line_phase,index_ls,re_phase) 
%
    pi_2 = pi*2.0;
    line_tmp_phase = line_phase;
%-----------
    diff_tmp = re_phase - line_phase;
%-----------
    index_bridge = zeros(1,length(line_phase)); 
    diff_bridge = zeros(1,length(line_phase));
    seg_mean_shift = zeros(1,length(index_ls));
    phi_g = 1.0*pi;
    kk = 1;
    for ii = 1:length(line_phase)
        if re_phase(ii) == 0;
            diff_tmp(ii) = 0;
        else
            index_bridge(kk) = ii;
            diff_bridge(kk) = diff_tmp(ii);
            kk = kk + 1;
        end
    end
    index_bridge(kk:length(line_phase)) = [];
    diff_bridge(kk:length(line_phase)) = [];
%------------------------------------------------------------------------
%   check the big long jump
    flag_display = 0;
    %index_p = find(diff_tmp > 1.2*pi);
    %index_n = find(diff_tmp < -1.2*pi);
    %
    %if isempty(index_pole_y)
    %else
    %   index_p
    %   index_n
    %   plot(line_phase,'o-g');hold on;
    %   plot(re_phase,'o-k');hold off;
    %   pause
    %end    
%------------------------------------------------------------------------
    %index_good = find(abs(diff_bridge)< pi)
    %pause
    %index_bridge(find(diff_bridge > 1.5*pi))
    %index_bridge(find(diff_bridge < -1.5*pi))
%------------------------------------------------------------------------
    seg_mean_shift = zeros(1,length(index_ls));
    for ii = 1:2:(length(index_ls) - 1)
    kk = 0;
    sum_diff = 0;
        for jj = index_ls(ii):index_ls(ii+1)
            if re_phase(jj) ~= 0
                sum_diff = sum_diff + diff_tmp(jj);
                kk = kk + 1;
            end
        end    
    %
        if kk > 0
            ave_diff = sum_diff/kk;
            seg_mean_shift(ii) = ave_diff;
            seg_mean_shift(ii+1) = ave_diff;           
        end
    end
%-------------------------------------------------------
    mean_trail = zeros(1,length(index_ls));
    for ii = 1:2:(length(index_ls) - 1)
        if index_ls(ii+1) - index_ls(ii) < 20;
            mm = 0;
            sum_mean = 0;
            for jj = index_ls(ii):index_ls(ii+1)
            sum_mean = sum_mean + line_phase(jj);
            mm = mm + 1;
            end 
            mean_trail(ii) = sum_mean/mm;
            mean_trail(ii + 1) = mean_trail(ii);
        else
            mm = 0;
            sum_mean = 0;
            for jj = index_ls(ii):(index_ls(ii)+10)
            sum_mean = sum_mean + line_phase(jj);
            mm = mm + 1;
            end 
            mean_trail(ii) = sum_mean/mm;
            %
            mm = 0;
            sum_mean = 0;
            for jj = (index_ls(ii+1) - 10):index_ls(ii + 1)
            sum_mean = sum_mean + line_phase(jj);
            mm = mm + 1;
            end 
            mean_trail(ii + 1) = sum_mean/mm;            
        end
    end
    diff_mean_trail = zeros(1, length(index_ls)/2);
    kk = 1;
    for ii = 2:2:(length(index_ls) - 2)
        diff_mean_trail(kk) = mean_trail(ii + 1) - mean_trail(ii);
        kk = kk + 1;
    end
    tmp_diff_mean_trail = find(abs(diff_mean_trail) > 1.0*pi);
%-------------------------------------------------------
  test_tmp = find(abs(seg_mean_shift) > phi_g);
  test_tmp_1 = find(abs(seg_mean_shift) == 0);
  flag_process = 1;
  if isempty(test_tmp) %|| isempty(tmp_diff_mean_trail)
     flag_process = 0;
  end
  %
  if flag_process == 1
    %
    for ii = 1:2:(length(index_ls) - 1)
       ave_diff = seg_mean_shift(ii);
       ls_s = index_ls(ii);
       ls_e = index_ls(ii+1);
       line_tmp_phase(ls_s:ls_e) = line_phase(ls_s:ls_e) + pi_2*round(ave_diff/(pi_2));  
    end
    %
    for ii = 2:2:(length(index_ls) - 2)
       in_s = index_ls(ii) + 1;
       in_e = index_ls(ii+1) -1;
       if in_e - in_s >= 0
           kk = 0; 
           sum_diff = 0;
            for jj = (index_ls(ii)-1):(index_ls(ii+1)+1);
                if re_phase(jj) ~= 0
                 sum_diff = sum_diff + diff_tmp(jj);
                 kk = kk + 1;
                end
            end
            if kk == 0
            ave_diff = (seg_mean_shift(ii) + seg_mean_shift(ii+1))/2;
            else
            ave_diff =  sum_diff/kk;%(seg_mean_shift(ii) + seg_mean_shift(ii+1))/2;
            end
            line_tmp_phase(in_s:in_e) = line_phase(in_s:in_e) + pi_2*round(ave_diff/(pi_2));
       end
    end
    %
    diff_check = re_phase - line_tmp_phase;
    index_check = find(abs(re_phase - line_tmp_phase) > pi_2*1000);
    for ii = 2:(length(index_check)-1)
        jj= index_check(ii);
        if re_phase(jj) ~= 0
            line_tmp_phase(jj) = line_tmp_phase(jj) + pi_2*round(diff_check(jj)/(pi_2*2));
        else
            if re_phase(jj-1) ~= 0 && re_phase(jj+1) ~= 0
                jump_mean = (diff_check(jj-1) + diff_check(jj) + diff_check(jj + 1))/3;
                line_tmp_phase(jj) = line_tmp_phase(jj) + pi_2*round(jump_mean/(pi_2*2));
            end
        end    
    end    
    %
        %figure(1)
        %plot(re_phase,'o-b');hold on;
        %plot(line_phase,'o-r');hold on;
        %index_ls
        %seg_mean_shift
        %mean_trail
        %plot(line_tmp_phase,'o-g');hold on;
        %plot((re_phase - line_phase),'o-k');hold off;
        %pause
  end%
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %     just for display check
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %flag_display = 0;
        if flag_display == 1
        figure(1)
        plot(re_phase,'o-b');hold on;
        plot(line_phase,'o-r');hold on;
        index_ls
        seg_mean_shift
        mean_trail
        plot(line_tmp_phase,'o-g');hold on;
        plot((re_phase - line_phase),'o-k');hold off;
        pause 
        end
%
%--------------------------------------------------------------------------