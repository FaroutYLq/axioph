function count_rate = auto_pulse_tagger(fine_iq_name, filenum, std_coeff, f_targ, power, CR_coeff,CR_thresh,make_plt, save_data)

    % Sample Pulse
    Ns = 2^9; % number of samples
    fs = 1E5;% sampling freq (in Hz)
    t0 =  1E-4; % where a test pulse starts (in seconds)
    tau = 3E-4; % time constant of pulse fall off (in seconds)
    [~, pulse]=PulseEg(Ns,fs,t0,tau);
    
    filebase= 'timeS';
    rates = zeros(size(filenum));

    for n = 1:numel(filenum)
        num = filenum(n);
        a = eConvertLVBinS( strcat(filebase, num2str(num), '.bin'));
        time = a(:,1);
        I=a(:,2);
        Q=a(:,3);
        theta = convert_IQ_to_theta_no_offset(I,Q,fine_iq_name);   

        % check if need to apply -1 to theta
        if abs(mean(theta)-min(theta)) > abs(mean(theta)-max(theta))
            theta = -theta;
        end
        
        signal = theta - mean(theta);
        avg_theta = movmean(signal,5);
    
        y = signal;
        
        yconv = conv(y,flip(pulse));
        yconv = yconv(Ns-1:end); % align
    
        signal = yconv;
    
        thresh = mean(signal) + std_coeff*std(signal);
        % second round if the data is noisy, and the baseline is too high. need to lower it in this case, by only using the data that is above the threshold to recalculate the threshold
        if 2*mean(abs(y)) < std(y)
            thresh_mask = signal < thresh;
            thresh = mean(signal) + std_coeff*std(signal(thresh_mask));
        end


        pulse_idx = find(signal<thresh);
        % minimum spacing between pulses is 20 samples; find the indicies of the first sample of each hit
        pulse_idx =  pulse_idx(diff(pulse_idx)>20); % cut down data size but not enought to miss a pulse
        
        pulse_idx2 = zeros(size(pulse_idx)); % index of the maximum
        CR_pulse_idx = zeros(size(pulse_idx)); % index of the CR
        for p = 1:size(pulse_idx,1)
            if pulse_idx(p)+61 < size(signal,1) % doesnt go off end of array
                % find the index of the minimum and maximum of the pulse
                [m,argmin] = min(signal(pulse_idx(p)+1:pulse_idx(p)+60));
                [M,argmax1] = max(signal(pulse_idx(p)+1:pulse_idx(p)+100)); % extend the range to find the maximum of the pulse
                [~,argmax2] = max(signal(pulse_idx(p)+1:pulse_idx(p)+60));

                if m > thresh && M < CR_coeff*std(signal)% whole hit is above threshold but not a CR
                    pulse_idx2(p) = pulse_idx(p) + argmax2 - 1; % set to max
                elseif m > thresh
                    CR_pulse_idx(p) = pulse_idx(p) + argmax1 - 1;
                end
            end
        end
        pulse_idx2 = nonzeros(pulse_idx2);
        CR_pulse_idx = nonzeros(CR_pulse_idx);
        pulse_idx =  pulse_idx2;    
    
       % align to top of peaks
        pulse_idx3 = zeros(size(pulse_idx));
        for p = 1:numel(pulse_idx)
            [~,argmax] = max(avg_theta( max([pulse_idx(p)-40 1]): min([pulse_idx(p)+40 5e6]) ));
            pulse_idx3(p) = pulse_idx(p) + argmax - 41;  
        end
        pulse_idx = pulse_idx3;
        mean_pulse_amp = mean(avg_theta(pulse_idx));
        std_pulse_amp = std(avg_theta(pulse_idx));

        pulse_idx3 = zeros(size(pulse_idx));
        pulse_idx3_CR = zeros(size(pulse_idx));
        for p = 1:numel(pulse_idx)
            [M,argmax] = max(avg_theta( max([pulse_idx(p)-40 1]): min([pulse_idx(p)+40 Ns]) ));
            % CR_coeff = 1;
            % CR_thresh = 1;
            if M >= min(CR_thresh, mean_pulse_amp + CR_coeff*std_pulse_amp)
                pulse_idx3_CR(p) = pulse_idx(p) + argmax - 41;
            else
                pulse_idx3(p) = pulse_idx(p) + argmax - 41;
            end
        end
        pulse_idx3 = nonzeros(pulse_idx3);
        CR_pulse_idx = cat(1,CR_pulse_idx,nonzeros(pulse_idx3_CR));
        pulse_idx = pulse_idx3;
        
         % extra layer of separating pulses from symmetric spikes
        fit_bool = 1;
        if fit_bool
            fit_arr = zeros(numel(pulse_idx),2); % each entry is [slope, y intercept]
            pulse_idx4 = zeros(size(pulse_idx));
            pulse_idx5 = zeros(size(pulse_idx));
            for p = 1:numel(pulse_idx)
                slice_size = 25;
                slice_before = avg_theta(pulse_idx(p)-slice_size:pulse_idx(p)); 
                for sl = 1:numel(slice_before)
                    if slice_before(sl) < 0
                        slice_before(sl) = 0;
                    end
                end
                linfit = polyfit(time(pulse_idx(p)-slice_size:pulse_idx(p)),slice_before,1);
                if linfit(1) > 75 % minimum slope
                    pulse_idx4(p) = pulse_idx(p);
                    fit_arr(p,1) = linfit(1);
                    fit_arr(p,2) = linfit(2);
                else
                    pulse_idx5(p) = pulse_idx(p);
                end
    
            end
            pulse_idx4 = nonzeros(pulse_idx4);
            pulse_idx5 = nonzeros(pulse_idx5);
            fit_arr( ~any(fit_arr,2), : ) = []; % remove zero rows
    
        end
    
        rates(n) = size(pulse_idx4,1)/100;
    
        if make_plt && n==1
            % pulse_idx = pulse_idx(avg_theta(pulse_idx) > 1.2e-6 / 8.597E-6);
            figure%('Position',[600 100 600 900]);
            hold on;
%             title(f_targ+"MHz, "+BB_T+"K," + stage_T + "mK")
            xlabel('time (s)')
            plot(time,avg_theta)
            scatter(time(pulse_idx),avg_theta(pulse_idx))
            scatter(time(CR_pulse_idx),avg_theta(CR_pulse_idx))


            for p = 1:numel(pulse_idx4)
                t_slice = time(pulse_idx4(p)-slice_size:pulse_idx4(p));
                plot(t_slice, fit_arr(p,1).*t_slice + fit_arr(p,2),"r")
            end
    
        end
    
        amplitudes{n} = avg_theta(pulse_idx4); 
        amplitudes_nonpulse{n} = avg_theta(pulse_idx5); 
        start_times{n} = time(pulse_idx4); % actually peak time
        start_times_nonpulse{n} = time(pulse_idx5); % actually peak time
        CR_times{n} = time(CR_pulse_idx);
    end
    
    if save_data
        save_path = "C:\Users\chalbert\Desktop\PRIMA_general\20250517_Cooldown\analysis\pulse_tagging\";
        % info_tag = f_targ+"MHz_"+BB_T+"K_"+num2str(stage_T)+"mK_"+std_coeff+"sigma";
        info_tag = power+"K_"+std_coeff+"sigma_"+f_targ+"Hz";
        save(save_path+"start_times_"+info_tag+".mat", 'start_times')
        save(save_path+"start_times_nonpulse_"+info_tag+".mat", 'start_times_nonpulse')
        save(save_path+"CR_times_"+info_tag+".mat", 'CR_times')
        save(save_path+"amplitudes_"+info_tag+".mat", 'amplitudes')
        save(save_path+"amplitudes_nonpulse_"+info_tag+".mat", 'amplitudes_nonpulse')
    end

    count_rate = mean(rates);
    % count_rate = rates;
    
end