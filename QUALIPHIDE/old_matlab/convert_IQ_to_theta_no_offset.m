function theta_TS = convert_IQ_to_theta_no_offset(I,Q,iq_fine_path)
    % Inputs:
    % I, Q - timestreams in I and Q
    % iq_fine_path - string path to fine scan of resonator
    %
    % Outputs:
    % theta_TS - timestream in angle with respect to IQ loop center. Should
    %               be used for relative theta shifts only because the
    %               theta=0 point is not chosen to be aligned to any
    %               particular feature of the loop in this script.

    finescan = load( iq_fine_path);
    
    [Icenter,Qcenter,IQ_rad,~] = circfit(finescan(:,2),finescan(:,3));

    lower_lim = atan2(finescan(1,3)-Qcenter,finescan(1,2)-Icenter);

    theta_TS = atan2(Q-Qcenter,I-Icenter); 

    jumps = theta_TS < lower_lim;
    theta_TS(jumps) = theta_TS(jumps) + 2*pi;

    

end
