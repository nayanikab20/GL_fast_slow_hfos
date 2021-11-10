% gui_findeventphase.m

function[event_phase]=gui_findeventphase(CORTEX,states,Mx,fn) 

%Band pass filter design:
    Wn1=[0.3/(fn/2) 300/(fn/2)]; 
    [b2,a2] = butter(3,Wn1); %0.3 to 300Hz
%Convert signal to 1 sec epochs.
        e_t=1;
        e_samples=e_t*(fn); %fs=1kHz
        ch=length(CORTEX);
        nc=floor(ch/e_samples); %Number of epochsw
        NC=[];
        for kk=1:nc
          NC(:,kk)= CORTEX(1+e_samples*(kk-1):e_samples*kk);
        end
        vec_bin=states;
        %Convert to 1 if NREM.
        vec_bin(vec_bin~=3)=0;
        vec_bin(vec_bin==3)=1;
        %Cluster one values:
        v2=ConsecutiveOnes(vec_bin);
        v_index=find(v2~=0);
        v_values=v2(v2~=0);
    for epoch_count=1:length(v_index)
    v{epoch_count,1}=reshape(NC(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);
    end
    V=cellfun(@(equis) filtfilt(b2,a2,equis), v ,'UniformOutput',false); %0.3 to 300Hz
    
% Bandpass filter design for Slow wave Oscillations
    Wn1=[0.5/(fn/2) 4/(fn/2)];
    [b1,a1]=butter(2, Wn1);
    V_SO=cellfun(@(equis) filtfilt(b1,a1,equis), v ,'UniformOutput',false); %0.5 to 4Hz
    
% Calculate phase of Slow wave oscillations
    V_SO_phase = cellfun(@(equis) mod(rad2deg(angle(hilbert(equis))),360) , V_SO,'UniformOutput',false); % Phases wrt SO
    
% Calculate phase of event wrt slow wave oscillations

event_timestamp=cellfun(@(equis1) uint16(equis1.*1000), Mx, 'UniformOutput', false); % transform event peak relative to epoch into timestamp
if ~isempty(Mx)
        event_phase =cellfun(@(equis1, equis2) equis1(equis2), V_SO_phase', event_timestamp , 'UniformOutput', false);
end

end



   
    