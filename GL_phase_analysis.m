% GL_phase_analysis.m
% Computes phase of the events with Slow wave oscillations (0.5-2 Hz)- HFO, spindles, Ripples 

% Requires 'load_me_first.mat' loaded first. 
% run first for rat 24 then 26 and 27

%% Find location
close all
dname=uigetdir([],'Select folder with Matlab data containing all rats.');
cd(dname)
%cd('/home/adrian/Documents/Plusmaze_downsampled')

%%
%Select rat number
opts.Resize = 'on';
opts.WindowStyle = 'modal';
opts.Interpreter = 'tex';
prompt=strcat('\bf Select a rat#. Options:','{ }',num2str(rats));
answer = inputdlg(prompt,'Input',[2 30],{''},opts);
Rat=str2num(answer{1});
cd(num2str(Rat))
tr=getfield(T,strcat('Rat',num2str(Rat)));%Thresholds 
%%
%Cortical regions.
yy={'PAR'};    
xx={'PFC'};  
%Sampling freq.
fn=1000;

%%
labelconditions=getfolder;
labelconditions=labelconditions.';
g=labelconditions;

multiplets=[{'singlets'} {'doublets'} {'triplets'} {'quatruplets'} {'pentuplets'} {'sextuplets'} {'septuplets'} {'octuplets'} {'nonuplets'}];
iii=1;

%% Select conditions and sessions
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
%xo
if sum(cell2mat(cellfun(@(equis1) contains(equis1,'nl'),g,'UniformOutput',false)))==1
g=g([find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{1}),g,'UniformOutput',false)))...
 find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{2}),g,'UniformOutput',false)))...
 find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{3}),g,'UniformOutput',false)))...
 find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{4}),g,'UniformOutput',false)))]);

else
%     error('Name issue')
end

if(Rat==24)
    PFC_spindle_phase_baseline=[]; PFC_spindle_phase_novelty=[]; PFC_spindle_phase_foraging=[]; PFC_spindle_phase_plusmaze =[]; 
    PPC_spindle_phase_baseline=[]; PPC_spindle_phase_novelty=[]; PPC_spindle_phase_foraging=[]; PPC_spindle_phase_plusmaze =[]; 
    slow_HFO_phase_baseline=[]; slow_HFO_phase_novelty=[]; slow_HFO_phase_foraging=[]; slow_HFO_phase_plusmaze =[]; 
    fast_HFO_phase_baseline=[]; fast_HFO_phase_novelty=[]; fast_HFO_phase_foraging=[]; fast_HFO_phase_plusmaze =[]; 
    hpc_ripple_phase_baseline=[];hpc_ripple_phase_novelty=[];hpc_ripple_phase_foraging=[];hpc_ripple_phase_plusmaze=[];
end
%%
f=waitbar(0,'Please wait...');
for k=1:length(g)
cd(g{k})

%% Find PFC events - Spindles
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

% PFC spindles
[spindle_pfc,~,~,Mx_pfc,timeasleep,sig_pfc,Ex_pfc,Sx_pfc]=gui_findspindlesZugaro(CORTEX,states,xx,multiplets,fn);

% phase of PFC spindles
PFC_spindle_phase=gui_findeventphase(CORTEX,states,Mx_pfc,fn);
PFC_spindle_phase_single = vertcat(PFC_spindle_phase{:});
eval(['PFC_spindle_phase_' labelconditions2{k} '=[PFC_spindle_phase_' labelconditions2{k} '; vertcat(PFC_spindle_phase{:}) ];']);

figure,
% Polar histogram
subplot(3,2,1)
polarhistogram(PFC_spindle_phase_single, 20)
heading=[labelconditions2{k} '-PFC Spindle'];
title(heading)

%% Find PPC events - Spindles, HFOS
par=dir(strcat('*',yy{1},'*.mat'));
if isempty(par)
    g=g(~contains(g,g{k}));
    cd ..
    progress_bar(k,length(g),f)
    break
end

par=par.name;
par=load(par);
par=getfield(par,yy{1});
par=par.*(0.195);

% PPC spindles
[spindle_ppc,~,~,Mx_ppc,timeasleep,sig_ppc,Ex_ppc,Sx_ppc]=gui_findspindlesZugaro(par,states,yy,multiplets,fn);

% phase of PPC spindles
PPC_spindle_phase=gui_findeventphase(par,states,Mx_ppc,fn);
PPC_spindle_phase_single = vertcat(PPC_spindle_phase{:});
eval(['PPC_spindle_phase_' labelconditions2{k} '=[PPC_spindle_phase_' labelconditions2{k} '; vertcat(PPC_spindle_phase{:}) ];']);

% Polar histogram
subplot(3,2,2)
polarhistogram(PPC_spindle_phase_single, 20)
heading=[labelconditions2{k} '-PPC Spindle'];
title(heading)


% Find PPC HFOs
[ripple_hfo,~,~,Mx_hfo,~,sig_hfo,Ex_hfo,Sx_hfo,~,~,~,~,~]=gui_findripples_Lisa(par,states,yy,tr,multiplets,fn);

%HFOs waveforms
si=sig_hfo(~cellfun('isempty',sig_hfo));
si=[si{:}];
% Group in slow and fast HFOs
[~,~,~,~,~,~,~,~,si_mixed,~]=hfo_specs(si,timeasleep,0,Rat,tr);
%% Separate bimodal distribution (Average freq)  in slow and fast HFOs
%g1: slow, g2:fast
%Initialize variables
Mx_hfo_g1=Mx_hfo;
Mx_hfo_g2=Mx_hfo;
Ex_hfo_g1=Ex_hfo;
Ex_hfo_g2=Ex_hfo;
Sx_hfo_g1=Sx_hfo;
Sx_hfo_g2=Sx_hfo;

row=si_mixed.i1;
cont=0;
for ll=1:length(Mx_hfo)

    if ~isempty(Mx_hfo{ll})

        for lll=1:length(Mx_hfo{ll})
            cont=cont+1;

            if ~ismember(cont,row)
                Mx_hfo_g1{ll}(lll)=NaN;
                Ex_hfo_g1{ll}(lll)=NaN;
                Sx_hfo_g1{ll}(lll)=NaN;
                
            else
                Mx_hfo_g2{ll}(lll)=NaN;
                Ex_hfo_g2{ll}(lll)=NaN;
                Sx_hfo_g2{ll}(lll)=NaN;

            end

        end
         Mx_hfo_g1{ll}=Mx_hfo_g1{ll}(~isnan(Mx_hfo_g1{ll}));
         Mx_hfo_g2{ll}=Mx_hfo_g2{ll}(~isnan(Mx_hfo_g2{ll}));

         Ex_hfo_g1{ll}=Ex_hfo_g1{ll}(~isnan(Ex_hfo_g1{ll}));
         Ex_hfo_g2{ll}=Ex_hfo_g2{ll}(~isnan(Ex_hfo_g2{ll}));
         Sx_hfo_g1{ll}=Sx_hfo_g1{ll}(~isnan(Sx_hfo_g1{ll}));
         Sx_hfo_g2{ll}=Sx_hfo_g2{ll}(~isnan(Sx_hfo_g2{ll}));
         
         
    end

end

% phase of PPC HFOs
% Slow
slow_HFO_phase=gui_findeventphase(par,states,Mx_hfo_g1,fn);
slow_HFO_phase_single = vertcat(slow_HFO_phase{:});
eval(['slow_HFO_phase_' labelconditions2{k} '=[slow_HFO_phase_' labelconditions2{k} '; vertcat(slow_HFO_phase{:}) ];']);

% Polar histogram
subplot(3,2,3)
polarhistogram(slow_HFO_phase_single, 20)
heading=[labelconditions2{k} '-Slow HFO'];
title(heading)

%fast
fast_HFO_phase=gui_findeventphase(par, states, Mx_hfo_g2,fn);
fast_HFO_phase_single = vertcat(fast_HFO_phase{:});
eval(['fast_HFO_phase_' labelconditions2{k} '=[fast_HFO_phase_' labelconditions2{k} '; vertcat(fast_HFO_phase{:}) ];']);

% Polar histogram
subplot(3,2,4)
polarhistogram(fast_HFO_phase_single, 20)
heading=[labelconditions2{k} '-Fast HFO'];
title(heading)

%% HPC Ripples
HPC=dir(strcat('*','HPC','*.mat'));
HPC=HPC.name;
HPC=load(HPC);
HPC=getfield(HPC,'HPC');
HPC=HPC.*(0.195);
%Find ripples
[ripple_nhpc,~,~,Mx_nhpc,~,sig_nhpc,Ex_nhpc,Sx_nhpc, ...
  ~,~,~,~,~,ripples_nhpc_vec, ti_cont, duration_nhpc_epoch...
  ]=gui_findripples_Lisa(HPC,states,{'HPC'},tr,multiplets,fn);

hpc_ripple_phase=gui_findeventphase(HPC, states, Mx_nhpc,fn);
hpc_ripple_phase_single = vertcat(hpc_ripple_phase{:});
eval(['hpc_ripple_phase_' labelconditions2{k} '=[hpc_ripple_phase_' labelconditions2{k} '; vertcat(hpc_ripple_phase{:}) ];']);

% Polar histogram
subplot(3,2,5)
polarhistogram(hpc_ripple_phase_single, 20)
heading=[labelconditions2{k} '-HPC Ripple'];
title(heading)


progress_bar(k,length(g),f)
cd ..    
end

xo
%% plotting polar histogram

%slow HFO
figure,
for k=1:4
    subplot(2,2,k)
    eval(['polarhistogram([slow_HFO_phase_' labelconditions2{k} '] , 20);'])
    heading=[labelconditions2{k} '-Slow HFO'];
    title(heading)
end

%fast HFO
figure,
for k=1:4
    subplot(2,2,k)
    eval(['polarhistogram([fast_HFO_phase_' labelconditions2{k} '] , 20);'])
    heading=[labelconditions2{k} '-fast HFO'];
    title(heading)
end

% PFC Spindles
figure,
for k=1:4
    subplot(2,2,k)
    eval(['polarhistogram([PFC_spindle_phase_' labelconditions2{k} '] , 20);'])
    heading=[labelconditions2{k} '-PFC spindle'];
    title(heading)
end

% PPC Spindles
figure,
for k=1:4
    subplot(2,2,k)
    eval(['polarhistogram([PPC_spindle_phase_' labelconditions2{k} '] , 20);'])
    heading=[labelconditions2{k} '-PPC spindle'];
    title(heading)
end

% HPC ripple
figure,
for k=1:4
    subplot(2,2,k)
    eval(['polarhistogram([hpc_ripple_phase_' labelconditions2{k} '] , 20);'])
    heading=[labelconditions2{k} '-HPC Ripple'];
    title(heading)
end