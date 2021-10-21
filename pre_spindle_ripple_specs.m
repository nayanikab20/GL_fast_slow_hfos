%GL_spindles_control
%Detects spindles and computes coocurrence with hfos and ripples.
%Controls by changing timestamps randomly to generate a null-distribution.
% Requires 'load_me_first.mat' loaded first. 

%% Find location
% close all
dname=uigetdir([],'Select folder with Matlab data containing all rats.');
cd(dname)
%cd('/home/adrian/Documents/Plusmaze_downsampled')

%%
%Select rat ID
opts.Resize = 'on';
opts.WindowStyle = 'modal';
opts.Interpreter = 'tex';
prompt=strcat('\bf Select a rat#. Options:','{ }',num2str(rats));
answer = inputdlg(prompt,'Input',[2 30],{''},opts);
Rat=str2num(answer{1});
cd(num2str(Rat))
%%
%Cortical regions.
yy={'PAR'};    
xx={'PFC'};  
%Sampling freq.          
fn=1000;
%% Get folder names
labelconditions=getfolder;
labelconditions=labelconditions.';
g=labelconditions;

multiplets=[{'singlets'} {'doublets'} {'triplets'} {'quatruplets'} {'pentuplets'} {'sextuplets'} {'septuplets'} {'octuplets'} {'nonuplets'}];
iii=1;

%%  Select conditions and sessions
    %Center figure.
    f=figure();
    movegui(gcf,'center');

    %Checkboxes
    Boxcheck = cell(1,4);
    for h1=1:length(labelconditions)
    boxcheck = uicontrol(f,'Style','checkbox','String',labelconditions{h1},'Position',[10 f.Position(4)-30*h1 400 20]);
    boxcheck.FontSize=11;
    boxcheck.Value=1;
    Boxcheck{h1}=boxcheck;   
    end

    set(f, 'NumberTitle', 'off', ...
        'Name', 'Select conditions');

    %Push button
    c = uicontrol;
    c.String = 'Continue';
    c.FontSize=10;
    c.Position=[f.Position(1)/3.5 c.Position(2)-10 f.Position(3)/2 c.Position(4)];

    %Callback
    c.Callback='uiresume(gcbf)';
    uiwait(gcf); 
    boxch=cellfun(@(x) get(x,'Value'),Boxcheck);
    clear Boxcheck
    close(f);
g={g{logical(boxch)}};    
if sum(cell2mat(cellfun(@(equis1) contains(equis1,'nl'),g,'UniformOutput',false)))==1
g=g([find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{1}),g,'UniformOutput',false)))...
 find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{2}),g,'UniformOutput',false)))...
 find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{3}),g,'UniformOutput',false)))...
 find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{4}),g,'UniformOutput',false)))]);

else
%    error('Name issue')
end

%Get thresholds for event detection.
tr=getfield(T,strcat('Rat',num2str(Rat)));%Thresholds 

ratno= ['rat', num2str(Rat)];

if exist('ts_hpc_pfc_pre')~=1
    
    ts_hpc_pfc_pre=struct();
    ts_pfc_hpc_pre=struct();
    ts_hpc_pfc_post=struct();
    ts_pfc_hpc_post=struct();
    
    duration=struct();
    amplitude=struct();

end


%%
f=waitbar(0,'Please wait...');
    for k=1:length(g)
        cd(g{k})
CORTEX=dir(strcat('*',xx{1},'*.mat'));
if isempty(CORTEX)
    g=g(~contains(g,g{k}));
    cd ..
    progress_bar(k,length(g),f)
    break
end
CORTEX=CORTEX.name;
CORTEX=load(CORTEX);
CORTEX=getfield(CORTEX,xx{1});
CORTEX=CORTEX.*(0.195);

A = dir('*states*.mat');
A={A.name};

if  ~isempty(A)
       cellfun(@load,A);
else
      error('No Scoring found')    
end

%% Find PFC Spindles
clear Ex_cortex Sx_cortex Mx_cortex
[ripple,~,~,Mx_cortex,timeasleep,sig_cortex,Ex_cortex,Sx_cortex, spindles_vec, ti_cont, ...
    start_epoch_cortex]=gui_findspindlesZugaro(CORTEX,states,xx,multiplets,fn);

si=sig_cortex(~cellfun('isempty',sig_cortex));
si=[si{:}];

%%
clear dur_cortex
% Find duration of spindles
dur_cortex= spindles_vec(:,3)-spindles_vec(:,1);
dur_cortex=dur_cortex.';
clear Sr_cortex Er_cortex

% Shuffling timestamps randomly 1000 times    
l=length(spindles_vec(:,1));

%PRE POST ANALYSIS
[Sx_pre,Mx_pre,Ex_pre,Sx_post,Mx_post,Ex_post]=cellfun(@(equis1,equis2,equis3) pre_post_spindle(equis1,equis2,equis3) ,Sx_cortex,Mx_cortex,Ex_cortex ,'UniformOutput',false);

% converting time stamp to global reference frame
Sx_pre_vec=cellfun(@(equis1, equis2) equis1+equis2, Sx_pre, mat2cell([0; start_epoch_cortex(1:end-1)]',1, ones(1,length(Sx_pre))), 'UniformOutput', false);
Sx_pre_vec=[Sx_pre_vec{:}];
Ex_pre_vec=cellfun(@(equis1, equis2) equis1+equis2, Ex_pre, mat2cell([0; start_epoch_cortex(1:end-1)]',1, ones(1,length(Sx_pre))), 'UniformOutput', false);
Ex_pre_vec=[Ex_pre_vec{:}];

Sx_post_vec=cellfun(@(equis1, equis2) equis1+equis2, Sx_post, mat2cell([0; start_epoch_cortex(1:end-1)]',1, ones(1,length(Sx_post))), 'UniformOutput', false);
Sx_post_vec=[Sx_post_vec{:}];
Ex_post_vec=cellfun(@(equis1, equis2) equis1+equis2, Ex_post, mat2cell([0; start_epoch_cortex(1:end-1)]',1, ones(1,length(Sx_post))), 'UniformOutput', false);
Ex_post_vec=[Ex_post_vec{:}];
  
%% HPC ripples

CORTEX=dir(strcat('*','HPC','*.mat'));
if isempty(CORTEX)
    g=g(~contains(g,g{k}));
    cd ..
    progress_bar(k,length(g),f)
    break
end
CORTEX=CORTEX.name;
CORTEX=load(CORTEX);
CORTEX=getfield(CORTEX,'HPC');
CORTEX=CORTEX.*(0.195);

A = dir('*states*.mat');
A={A.name};
if  ~isempty(A)
       cellfun(@load,A);
else
      error('No Scoring found')    
end
% Find ripples
[ripple_nhpc,~,~,Mx_nhpc,~,sig_nhpc,Ex_nhpc,Sx_nhpc, ...
  ~,~,~,~,~,ripples_nhpc_vec, ti_cont, duration_nhpc_epoch...
  ]=gui_findripples_Lisa(CORTEX,states,{'HPC'},tr,multiplets,fn);

%% Coocur PFC spindle and HPC ripples

[ts_hpc_pfc.(ratno).(labelconditions2{k}),ts_pfc_hpc.(ratno).(labelconditions2{k}),count_hpc_pfc(k),count_pfc_hpc(k)]=cooccurrence_vec(ripples_nhpc_vec(:,1), ripples_nhpc_vec(:,3), spindles_vec(:,1), spindles_vec(:,3) );
[ts_hpc_pfc_pre.(ratno).(labelconditions2{k}),ts_pfc_hpc_pre.(ratno).(labelconditions2{k}),count_hpc_pfc_pre(k),count_pfc_pre_hpc(k)]=cooccurrence_vec(ripples_nhpc_vec(:,1), ripples_nhpc_vec(:,3), Sx_pre_vec',Ex_pre_vec' );
[ts_hpc_pfc_post.(ratno).(labelconditions2{k}),ts_pfc_hpc_post.(ratno).(labelconditions2{k}),count_hpc_pfc_post(k),count_pfc_post_hpc(k)]=cooccurrence_vec(ripples_nhpc_vec(:,1), ripples_nhpc_vec(:,3), Sx_post_vec',Ex_post_vec' );

%% specs of the ripples
%duration

duration.(ratno).([labelconditions2{k}, '_pre']) = ripples_nhpc_vec(ts_hpc_pfc_pre.(ratno).(labelconditions2{k}), 3)-ripples_nhpc_vec(ts_hpc_pfc_pre.(ratno).(labelconditions2{k}), 1);
duration.(ratno).([labelconditions2{k}, '_pre_mean'])=mean(duration.(ratno).([labelconditions2{k}, '_pre']));
duration.(ratno).(labelconditions2{k}) = ripples_nhpc_vec(ts_hpc_pfc.(ratno).(labelconditions2{k}), 3)-ripples_nhpc_vec(ts_hpc_pfc.(ratno).(labelconditions2{k}), 1);
duration.(ratno).([labelconditions2{k}, '_mean'])=mean(duration.(ratno).(labelconditions2{k}));


% amplitude
amplitude.(ratno).([labelconditions2{k}, '_pre']) = ripples_nhpc_vec(ts_hpc_pfc_pre.(ratno).(labelconditions2{k}), 2);
amplitude.(ratno).([labelconditions2{k}, '_pre_mean'])=mean(amplitude.(ratno).([labelconditions2{k}, '_pre']));
amplitude.(ratno).(labelconditions2{k}) = ripples_nhpc_vec(ts_hpc_pfc.(ratno).(labelconditions2{k}), 2);
amplitude.(ratno).([labelconditions2{k}, '_mean'])=mean(amplitude.(ratno).(labelconditions2{k}));

progress_bar(k,length(g),f)
    cd ..    
end
    
xo
for rat=1:3

    ratno=['rat' num2str(rats(rat))]
    
    figure;
    for k=1:4    
    % histograms
    edges=0:0.01:0.6;
    subplot(2,4,k)
    histogram(duration.(ratno).([labelconditions2{k}, '_pre']), edges, 'normalization', 'probability')
    hold on;
    histogram(duration.(ratno).(labelconditions2{k}), edges, 'normalization', 'probability')
    hold off;
    title(labelconditions2{k})
    xlabel('duration')
    ylabel('# ripples')
    legend('pre spindle ripple', 'spindle ripple')
    
    edges2=0:500:8000;
    subplot(2,4,4+k)
    histogram(amplitude.(ratno).([labelconditions2{k}, '_pre']), edges2, 'normalization', 'probability')
    hold on;
    histogram(amplitude.(ratno).(labelconditions2{k}), edges2, 'normalization', 'probability')
    hold off;
    title(labelconditions2{k})
    xlabel('amplitude')
    ylabel('# ripples')
    legend('pre spindle ripple', 'spindle ripple')
    end
end