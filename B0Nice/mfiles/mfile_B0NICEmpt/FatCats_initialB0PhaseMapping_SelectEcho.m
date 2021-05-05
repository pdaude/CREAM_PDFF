%
function [algoParams] = FatCats_initialB0PhaseMapping_SelectEcho(algoParams)
%
    pi_2 = 2*pi;
    TE_seq = algoParams.TE_seq;
    phase_echo = algoParams.Phase4Echo_ori;
%
if algoParams.echo_selection == 1
    close all
    plot(algoParams.Phase4Echo_ori,'r--*');
    xlabel('index of echo') % x-axis label
    ylabel('phase shift between fat- and water') % y-axis label
            flag_input = 0;
            while flag_input == 0
                clear test_prompt;
                disp('flag_B0 = 1, initial B0 map will be calculated from the HP between two selected echoes');
                disp('flag_B0 = 2, only for case #17');
                %disp('flag_B0 = 3, dualecho acqusition');
                test_prompt = input('please input the value of flag_B0 (mostlikely == 1) : ');
                if ~isempty(test_prompt)
                flag_B0 = test_prompt;
                else
                flag_B0 = 1;    
                end
                
                test_prompt = input('please input the index of echo at the early side (index_start) : ');
                if ~isempty(test_prompt)
                index_start = test_prompt;
                end
                
                test_prompt = input('please input the index of echo at the later side (index_end) : ');
                if ~isempty(test_prompt)
                index_end = test_prompt;
                end
                
                flag_input = 1;
            end
    index_B0 = [flag_B0 index_start index_end];
else

    diff_echo = diff(unwrap(phase_echo));
    diff_echo = sum(abs(diff_echo))/length(diff_echo);
    ratio_2piOverDiff = round(abs(pi_2/abs(diff_echo)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flag_B0 = 0;index_start = 0;index_end = 0;
    if length(phase_echo) >= 3
    %
        if (1+ratio_2piOverDiff) == length(TE_seq)
        flag_B0 = 1;
        index_start = 1;
        index_end = index_start + ratio_2piOverDiff; 
        elseif (1+ratio_2piOverDiff) < length(TE_seq)
            flag_B0 = 1;
            if abs(phase_echo(1)) < 0.9*pi && abs(phase_echo(1 + ratio_2piOverDiff)) < 0.9*pi
            index_start = 1;
            index_end = index_start + ratio_2piOverDiff;
            else
            index_start = 2;
            index_end = index_start + ratio_2piOverDiff;        
            end
        elseif ratio_2piOverDiff == 3 && length(TE_seq) == 3
        flag_B0 = 2;
        index_start = 1;
        index_end = length(TE_seq);            
        else
        flag_B0 = 1;
        index_start = 1;
        index_end = length(TE_seq);
        end
    %plot(phase_echo./pi,'o-r');hold on;
    %plot(unwrap(phase_echo)./pi,'o-k');hold off;
    end
%
    index_B0 = [flag_B0 index_start index_end];
    %end
    if index_end < length(TE_seq)
        if abs(phase_echo(index_end+1) - phase_echo(index_start)) < abs(phase_echo(index_end) - phase_echo(index_start))
        index_end = index_end + 1;
        index_B0 = [flag_B0 index_start index_end];
        end
    end
end
%
index_B0

if index_B0(1) == 1
    delta_TE4B0 = TE_seq(index_B0(3)) - TE_seq(index_B0(2)); 
    delta_PhiIni = (phase_echo(index_B0(3)) - phase_echo(index_B0(2)));
    delta_PhiIni = delta_PhiIni - pi_2*round(delta_PhiIni/pi_2);
elseif index_B0(1) == 2 %for case 17
    delta_TE4B0 = (TE_seq(index_B0(3))+TE_seq(index_B0(2)+1) - 2*TE_seq(index_B0(2)));
     delta_PhiIni = phase_echo(index_B0(3)) + phase_echo(index_B0(3)-1)- 2*phase_echo(index_B0(2));    
    if abs(delta_PhiIni) > pi
    delta_PhiIni = delta_PhiIni - sign(delta_PhiIni)*pi_2;       
    end   
elseif index_B0(1) == 3 % for dual echo acquistion
    delta_TE4B0 = (TE_seq(index_B0(3)) - TE_seq(index_B0(2)))*2;
    delta_PhiIni = 2*(phase_echo(index_B0(3)) - phase_echo(index_B0(2)));
    delta_PhiIni = delta_PhiIni - pi_2*round(delta_PhiIni/pi_2);    
end
%
algoParams.index_B0 = index_B0;
algoParams.delta_TE4B0 = delta_TE4B0;
algoParams.error_bias = delta_PhiIni/(2*pi*algoParams.delta_TE4B0);    
%
phase_echo = zeros(1,length(TE_seq));
mag_echo = zeros(1,length(TE_seq));
model_f = algoParams.model_f - algoParams.error_bias;
    for index_echo = 1:length(TE_seq)
    TE_tmp = TE_seq(index_echo);
    model_complex = algoParams.model_r.*exp(1i.*(2*pi*TE_tmp.*model_f));
    model_complex_sum = sum(model_complex);
    phase_echo(index_echo) = angle(model_complex_sum);
    mag_echo(index_echo) = abs(model_complex_sum);
    end
%
algoParams.phase_echo = phase_echo;
%

%