function [ripple,RipFreq,rip_duration,Mx,timeasleep,sig,Ex,Sx,ripple_multiplets,RipFreq_multiplets,rip_duration_multiplets,sig_multiplets,M_multiplets,V,Mono]=gui_findripples(CORTEX,states,xx,tr,multiplets,fn)
    %Band pass filter design:
    Wn1=[100/(fn/2) 300/(fn/2)]; % Cutoff=100-300 Hz
    [b1,a1] = butter(3,Wn1,'bandpass'); %Filter coefficients
    Wn1=[320/(fn/2)]; % Cutoff=320 Hz
    [b2,a2] = butter(3,Wn1); %Filter coefficients
%Rearrange signal in 1 sec epochs.
        e_t=1;
        e_samples=e_t*(fn); %fs=1kHz
        ch=length(CORTEX);
        nc=floor(ch/e_samples); %Number of epochsw
        NC=[];
        for kk=1:nc
          NC(:,kk)= CORTEX(1+e_samples*(kk-1):e_samples*kk);
        end
        %Find if epoch is NREM (state=3)
        vec_bin=states;
        vec_bin(vec_bin~=3)=0;
        vec_bin(vec_bin==3)=1;
        %Cluster one values:
        v2=ConsecutiveOnes(vec_bin);
        v_index=find(v2~=0);
        v_values=v2(v2~=0);
    %Extract NREM epochs    
    for epoch_count=1:length(v_index)
    v{epoch_count,1}=reshape(NC(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);
    end
    %Filter
    V=cellfun(@(equis) filtfilt(b2,a2,equis), v ,'UniformOutput',false);
    Mono=cellfun(@(equis) filtfilt(b1,a1,equis), V ,'UniformOutput',false);
    %Total amount of NREM time:
    timeasleep=sum(cellfun('length',V))*(1/fn)/60; % In minutes
    %NREM epochs without OpenEphys BitVolt factor
    signal2=cellfun(@(equis) times((1/0.195), equis)  ,Mono,'UniformOutput',false);
    %Timestamps for ripple/hfo detector
    ti=cellfun(@(equis) reshape(linspace(0, length(equis)-1,length(equis))*(1/fn),[],1) ,signal2,'UniformOutput',false);
    
    
    %% Zugaro detection algorithm
    
        %If HPC find ripples, else find hfos.
    %Finds start, end and peak timestamp.
    if strcmp(xx{1},'HPC')
        % Finding number of ripples in the dataset
        Concat_signal=vertcat(signal2{:});
        ti_cont=(1:length(Concat_signal))./1000;    
        Concat_input = [ti_cont' Concat_signal];
        [ripples, ~,~] =FindRipples_Zugaro(Concat_input, 'durations', [50, 30, 100], 'thresholds', [1 4] );

        %finding number of ripples per epoch

        %duration of each epoch of non-rem sleep
        duration_epoch=cellfun(@length, ti)/1000; 
        duration_epoch_cumsum=cumsum(duration_epoch);

        for f=1:length(duration_epoch_cumsum)

            if(f==1)
                vec=find(ripples(:,1)>=0 & ripples(:,3)<=duration_epoch_cumsum(f));
            elseif(f>1)
                vec=find(ripples(:,1)>=duration_epoch_cumsum(f-1) & ripples(:,3)<=duration_epoch_cumsum(f));
            end

            ripples_per_epoch(f)=length(vec);
            % assigning each spindle an epoch number
            ripples(vec,5)=f;
            %assigning each spindle a local start and stop time wrt the epoch
            if(f==1)
                ripples(vec, 6:8)=ripples(vec,1:3);
            elseif(f>1)
                ripples(vec, 6:8)=ripples(vec,1:3)-duration_epoch_cumsum(f-1);
            end

            Sx(f)={(ripples(vec,6))'};
            Ex(f)={(ripples(vec,8))'};
            Mx(f)= {(ripples(vec,7))'};
        end
    
    else
        [Sx,Ex,Mx] =cellfun(@(equis1,equis2) findHFOs(equis1, equis2, tr(2), (tr(2))*(1/2), [] ), signal2,ti,'UniformOutput',false);
    end
       
    for l=1:length(Mx)
         hfo_sequence=ConsecutiveOnes(diff(Mx{l})<=0.300);

         for ll=1:length(multiplets)
             eval([multiplets{ll} '=(hfo_sequence=='  num2str(ll-1) ');'])
             cont=1;
             M_multiplets.(multiplets{ll}){l}=[];
             S_multiplets.(multiplets{ll}){l}=[];
             E_multiplets.(multiplets{ll}){l}=[];
             
             while cont<=ll
                 eval(['Sx_' multiplets{ll} '_' num2str(cont) '{l}=Sx{l}(find(' multiplets{ll} ')+(cont-1));'])
                 eval(['Mx_' multiplets{ll} '_' num2str(cont) '{l}=Mx{l}(find(' multiplets{ll} ')+(cont-1));'])
                 eval(['Ex_' multiplets{ll} '_' num2str(cont) '{l}=Ex{l}(find(' multiplets{ll} ')+(cont-1));'])
                 Mx_multiplets.(multiplets{ll}).(strcat('m_',num2str(cont))){l}=eval(['Mx_' multiplets{ll} '_' num2str(cont) '{l}']);
                 M_multiplets.(multiplets{ll}){l}=eval(['sort([M_multiplets.(multiplets{ll}){l} ' ' Mx_' multiplets{ll} '_' num2str(cont) '{l}])']); % Combined consecutive multiplets    

                 Sx_multiplets.(multiplets{ll}).(strcat('m_',num2str(cont))){l}=eval(['Sx_' multiplets{ll} '_' num2str(cont) '{l}']);
                 S_multiplets.(multiplets{ll}){l}=eval(['sort([S_multiplets.(multiplets{ll}){l} ' ' Sx_' multiplets{ll} '_' num2str(cont) '{l}])']); % Combined consecutive multiplets    
                 
                 Ex_multiplets.(multiplets{ll}).(strcat('m_',num2str(cont))){l}=eval(['Ex_' multiplets{ll} '_' num2str(cont) '{l}']);
                 E_multiplets.(multiplets{ll}){l}=eval(['sort([E_multiplets.(multiplets{ll}){l} ' ' Ex_' multiplets{ll} '_' num2str(cont) '{l}])']); % Combined consecutive multiplets    

                 
                 %                   eval([  'clear' ' ' 'Mx_' multiplets{ll} '_' num2str(cont)])
                 cont=cont+1;
             end
         end
         
         
            %Correct singlets values
            multiplets_stamps=[];
            S_multiplets_stamps=[];
            E_multiplets_stamps=[];

            for i=2:length(multiplets)
                %multiplets_stamps=[multiplets_stamps ([M_multiplets.(multiplets{i}){:}])]; %[M_multiplets.doublets{:}]
                multiplets_stamps=[multiplets_stamps (M_multiplets.(multiplets{i}){l})];
                S_multiplets_stamps=[S_multiplets_stamps (S_multiplets.(multiplets{i}){l})];
                E_multiplets_stamps=[E_multiplets_stamps (E_multiplets.(multiplets{i}){l})];

            end

            M_multiplets.singlets(l)={ Mx{l}(~ismember(Mx{l}, multiplets_stamps))};
            Mx_multiplets.singlets.m_1(l)={Mx{l}(~ismember(Mx{l}, multiplets_stamps))};

            S_multiplets.singlets(l)={ Sx{l}(~ismember(Sx{l}, S_multiplets_stamps))};
            Sx_multiplets.singlets.m_1(l)={Sx{l}(~ismember(Sx{l}, S_multiplets_stamps))};

            E_multiplets.singlets(l)={ Ex{l}(~ismember(Ex{l}, E_multiplets_stamps))};
            Ex_multiplets.singlets.m_1(l)={Ex{l}(~ismember(Ex{l}, E_multiplets_stamps))};
            
    end
    
%     for l=1:length(Mx)
%          hfo_sequence=ConsecutiveOnes(diff(Mx{l})<=0.300);
% 
%          for ll=1:length(multiplets)
%              eval([multiplets{ll} '=(hfo_sequence=='  num2str(ll-1) ');'])
%              eval(['Sx_' multiplets{ll} '_1{l}=Sx{l}(find(' multiplets{ll} '));'])
%              eval(['Ex_' multiplets{ll} '_1{l}=Ex{l}(find(' multiplets{ll} '));'])
% 
%          end
%     end

% Get traces of events detected    
    for l=1:length(Sx)
         sig{l}=getsignal(Sx,Ex,ti,Mono,l);
%        sig{l}=getsignal(Sx,Ex,ti,V,l);
        for ll=1:length(multiplets)
            eval(['sig_' multiplets{ll} '_1{l}=getsignal(Sx_multiplets.' multiplets{ll} '.m_1, Ex_multiplets.' multiplets{ll} '.m_' num2str(ll) ',ti,Mono,l);'])
%            eval(['sig_multiplets.(' 'multiplets{ll}' ')=getsignal(Sx_' multiplets{ll} '_1,Ex_' multiplets{ll} '_1,ti,Mono,l);'])
                        eval(strcat('sig_',multiplets{ll},'_1=sig_',multiplets{ll},'_1.'';'))
              eval(['sig_multiplets.(' 'multiplets{ll}' ')=sig_' multiplets{ll},'_1.'';'])  
              eval(['clear' ' ' 'sig_' multiplets{ll} '_1'])
        end
    end

    sig=sig.';

    [ripple, RipFreq,rip_duration]=hfo_count_freq_duration(Sx,Ex,timeasleep);
    
    for ll=1:length(multiplets)
        eval(['[ripple_multiplets.(' 'multiplets{ll}' '), RipFreq_multiplets.(' 'multiplets{ll}' '),rip_duration_multiplets.(' 'multiplets{ll}' ')]=hfo_count_freq_duration(Sx_multiplets.' multiplets{ll} '.m_1, Ex_multiplets.' multiplets{ll} '.m_' num2str(ll) ',timeasleep);']);
    end

end