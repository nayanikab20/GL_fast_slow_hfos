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

tic
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

sec=1
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

Sr_cortex=zeros(1000,l);
parfor r=1:1000
    
     ti_rand=ti_cont(randperm(size(ti_cont,2)));
     Sr_cortex(r,:)=ti_rand(find(ismember(cast(ti_cont*1000, 'uint32'), cast(spindles_vec(:,1)*1000, 'uint32') )));
     Er_cortex(r,:)=Sr_cortex(r,:)+ dur_cortex;
end

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
      
sec =2
%% High Frequency oscillations

CORTEX=dir(strcat('*','PAR','*.mat'));
if isempty(CORTEX)
    g=g(~contains(g,g{k}));
    cd ..
    progress_bar(k,length(g),f)
    break
end
CORTEX=CORTEX.name;
CORTEX=load(CORTEX);
CORTEX=getfield(CORTEX,'PAR');
CORTEX=CORTEX.*(0.195);

A = dir('*states*.mat');
A={A.name};
if  ~isempty(A)
       cellfun(@load,A);
else
      error('No Scoring found')    
end
% Detect HFOs
[ripple_sphfo,~,~,Mx_cortex_sphfo,~,sig_cortex_sphfo,Ex_cortex_sphfo,Sx_cortex_sphfo,...
 ~,~,~,~,~,ripple_sphfo_vec ,ti_sphfo_cont, duration_sphfo_epoch...
  ]=gui_findripples_Lisa(CORTEX,states,{'PAR'},tr,multiplets,fn);

si=sig_cortex_sphfo(~cellfun('isempty',sig_cortex_sphfo));
si=[si{:}];

%Group events in slow and fast HFOs.
[~,~,~,~,~,~,~,~,si_mixed,~]=hfo_specs(si,timeasleep,0,Rat,tr);

%% Determining slow and fast HFOs timestamps.
%Initializing variables
Mx_cortex_g1=Mx_cortex_sphfo;
Mx_cortex_g2=Mx_cortex_sphfo;
Ex_cortex_g1=Ex_cortex_sphfo;
Ex_cortex_g2=Ex_cortex_sphfo;
Sx_cortex_g1=Sx_cortex_sphfo;
Sx_cortex_g2=Sx_cortex_sphfo;

row=si_mixed.i1;
cont=0;
for ll=1:length(Mx_cortex_sphfo)

    if ~isempty(Mx_cortex_sphfo{ll})

        for lll=1:length(Mx_cortex_sphfo{ll})
            cont=cont+1;
 
            if ~ismember(cont,row)
                Mx_cortex_g1{ll}(lll)=NaN;
                Ex_cortex_g1{ll}(lll)=NaN;
                Sx_cortex_g1{ll}(lll)=NaN;
                
            else
                Mx_cortex_g2{ll}(lll)=NaN;
                Ex_cortex_g2{ll}(lll)=NaN;
                Sx_cortex_g2{ll}(lll)=NaN;

            end

        end
         Mx_cortex_g1{ll}=Mx_cortex_g1{ll}(~isnan(Mx_cortex_g1{ll}));
         Mx_cortex_g2{ll}=Mx_cortex_g2{ll}(~isnan(Mx_cortex_g2{ll}));

         Ex_cortex_g1{ll}=Ex_cortex_g1{ll}(~isnan(Ex_cortex_g1{ll}));
         Ex_cortex_g2{ll}=Ex_cortex_g2{ll}(~isnan(Ex_cortex_g2{ll}));
         Sx_cortex_g1{ll}=Sx_cortex_g1{ll}(~isnan(Sx_cortex_g1{ll}));
         Sx_cortex_g2{ll}=Sx_cortex_g2{ll}(~isnan(Sx_cortex_g2{ll}));
         
         
    end

end

% converting cell wise local slow and fast hfo st, m, sp to global vectorised st, m, sp.
duration_epoch_cumsum=cumsum(duration_sphfo_epoch);
        duration_epoch_cumsum=[0; duration_epoch_cumsum];
        Sx_sphfo_g1=[]; Mx_sphfo_g1=[]; Ex_sphfo_g1=[]; Sx_sphfo_g2=[]; Mx_sphfo_g2=[]; Ex_sphfo_g2=[];
        for m=1:length(duration_epoch_cumsum)-1
            if ~isempty(Sx_cortex_g1(m))
                Sx_sphfo_g1= [Sx_sphfo_g1 cell2mat(Sx_cortex_g1(m))+ duration_epoch_cumsum(m)];
                Mx_sphfo_g1= [Mx_sphfo_g1 cell2mat(Mx_cortex_g1(m))+ duration_epoch_cumsum(m)];
                Ex_sphfo_g1= [Ex_sphfo_g1 cell2mat(Ex_cortex_g1(m))+ duration_epoch_cumsum(m)];
               
            end
             if ~isempty(Sx_cortex_g2(m))
                Sx_sphfo_g2= [Sx_sphfo_g2 cell2mat(Sx_cortex_g2(m))+ duration_epoch_cumsum(m)];
                Mx_sphfo_g2= [Mx_sphfo_g2 cell2mat(Mx_cortex_g2(m))+ duration_epoch_cumsum(m)];
                Ex_sphfo_g2= [Ex_sphfo_g2 cell2mat(Ex_cortex_g2(m))+ duration_epoch_cumsum(m)];
               
            end
        end

%% Coocur PFC spindle and hfos
[~,~, count_cohfos_pfc_g1(k), count_pfc_cohfos_g1(k)]= cooccurrence_vec( Sx_sphfo_g1,Ex_sphfo_g1, spindles_vec(:,1), spindles_vec(:,3));
[~,~, count_cohfos_pfc_g2(k), count_pfc_cohfos_g2(k)] = cooccurrence_vec( Sx_sphfo_g2,Sx_sphfo_g2, spindles_vec(:,1), spindles_vec(:,3) ) ;
[~,~, count_cohfos_pfc_pre_g1(k), count_pfc_pre_cohfos_g1(k)]= cooccurrence_vec( Sx_sphfo_g1,Ex_sphfo_g1, Sx_pre_vec',Ex_pre_vec');
[~,~, count_cohfos_pfc_pre_g2(k), count_pfc_pre_cohfos_g2(k)] = cooccurrence_vec( Sx_sphfo_g2,Sx_sphfo_g2, Sx_pre_vec',Ex_pre_vec' ) ;
[~,~, count_cohfos_pfc_post_g1(k), count_pfc_post_cohfos_g1(k)]= cooccurrence_vec( Sx_sphfo_g1,Ex_sphfo_g1, Sx_post_vec',Ex_post_vec' );
[~,~, count_cohfos_pfc_post_g2(k), count_pfc_post_cohfos_g2(k)] = cooccurrence_vec( Sx_sphfo_g2,Sx_sphfo_g2, Sx_post_vec', Ex_post_vec' ) ;

%% Coocur shuffled PFC spindle and original hfos.

parfor r=1:1000
    
    [~,~, count_cohfos_pfc_rand_g1(k,r), count_pfc_cohfos_rand_g1(k,r)]= cooccurrence_vec( Sx_sphfo_g1, Ex_sphfo_g1, Sr_cortex(r,:)', Er_cortex(r,:)' );
    [~,~, count_cohfos_pfc_rand_g2(k,r), count_pfc_cohfos_rand_g2(k,r)]= cooccurrence_vec( Sx_sphfo_g2, Ex_sphfo_g2, Sr_cortex(r,:)', Er_cortex(r,:)' );
    
end


sec=3
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

[~,~,count_hpc_pfc(k),count_pfc_hpc(k)]=cooccurrence_vec(ripples_nhpc_vec(:,1), ripples_nhpc_vec(:,3), spindles_vec(:,1), spindles_vec(:,3) );
[~,~,count_hpc_pfc_pre(k),count_pfc_pre_hpc(k)]=cooccurrence_vec(ripples_nhpc_vec(:,1), ripples_nhpc_vec(:,3), Sx_pre_vec',Ex_pre_vec' );
[~,~,count_hpc_pfc_post(k),count_pfc_post_hpc(k)]=cooccurrence_vec(ripples_nhpc_vec(:,1), ripples_nhpc_vec(:,3), Sx_post_vec',Ex_post_vec' );

%% Coocur shuffled PFC spindle and HPC ripples.

%Shuffle ripples timestamps
 parfor r=1:1000
     
     [~,~, count_hpc_pfc_rand(k,r), count_pfc_hpc_rand(k,r)]=cooccurrence_vec(ripples_nhpc_vec(:,1), ripples_nhpc_vec(:,3), Sr_cortex(r,:)', Er_cortex(r,:)' ) ;
 end


% %% COOCUR PFC SPINDLE-HPC MULTIPLETS
% [out_pfc]=coccur_multiplets(cohfos2_hpc);
% Out_PFC{k}=out_pfc;

sec=4
%% PARIETAL SPINDLES.

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
% Find Parietal spindle
[ripple,RipFreq,rip_duration,Mx_par,timeasleep,sig_par,Ex_par,Sx_par, spindles_ppc_vec, ti_ppc_cont, start_epoch_ppc]=gui_findspindlesZugaro(par,states,yy,multiplets,fn);


si=sig_par(~cellfun('isempty',sig_par));
si=[si{:}];

[x,y,z,~,~,~,l,p,si_mixed,th]=hfo_specs(si,timeasleep,0,Rat,tr);
%PRE POST ANALYSIS
[Sx_pre,Mx_pre,Ex_pre,Sx_post,Mx_post,Ex_post]=cellfun(@(equis1,equis2,equis3) pre_post_spindle(equis1,equis2,equis3) ,Sx_par,Mx_par,Ex_par ,'UniformOutput',false);

% converting time stamp to global reference frame
Sx_pre_vec=cellfun(@(equis1, equis2) equis1+equis2, Sx_pre, mat2cell([0; start_epoch_ppc(1:end-1)]',1, ones(1,length(Sx_pre))), 'UniformOutput', false);
Sx_pre_vec=[Sx_pre_vec{:}];
Ex_pre_vec=cellfun(@(equis1, equis2) equis1+equis2, Ex_pre, mat2cell([0; start_epoch_ppc(1:end-1)]',1, ones(1,length(Sx_pre))), 'UniformOutput', false);
Ex_pre_vec=[Ex_pre_vec{:}];

Sx_post_vec=cellfun(@(equis1, equis2) equis1+equis2, Sx_post, mat2cell([0; start_epoch_ppc(1:end-1)]',1, ones(1,length(Sx_post))), 'UniformOutput', false);
Sx_post_vec=[Sx_post_vec{:}];
Ex_post_vec=cellfun(@(equis1, equis2) equis1+equis2, Ex_post, mat2cell([0; start_epoch_ppc(1:end-1)]',1, ones(1,length(Sx_post))), 'UniformOutput', false);
Ex_post_vec=[Ex_post_vec{:}]; 
%%
%Find duration of spindles
clear dur_hpc dur_par
dur_par= spindles_ppc_vec(:,3)-spindles_ppc_vec(:,1);
dur_par=dur_par.';
clear Sr_par Er_par

% Shuffling timestamps of PAR Spindles randomly 1000 times    
l=length(spindles_ppc_vec(:,1));

Sr_par=zeros(1000,l);
parfor r=1:1000
    
     ti_rand=ti_ppc_cont(randperm(size(ti_ppc_cont,2)));
     Sr_par(r,:)=ti_rand(find(ismember(cast(ti_ppc_cont*1000, 'uint32'), cast(spindles_ppc_vec(:,1)*1000, 'uint32') )));
     Er_par(r,:)=Sr_par(r,:)+ dur_par;
end
   

%% Find coocurrent PAR spindles and hfos

[~,~, count_cohfos_par_g1(k), count_par_cohfos_g1(k)] = cooccurrence_vec( Sx_sphfo_g1, Ex_sphfo_g1, spindles_ppc_vec(:,1), spindles_ppc_vec(:,3)) ;
[~, ~, count_cohfos_par_g2(k), count_par_cohfos_g2(k)]= cooccurrence_vec( Sx_sphfo_g2, Ex_sphfo_g2, spindles_ppc_vec(:,1), spindles_ppc_vec(:,3) );
[~,~, count_cohfos_par_pre_g1(k), count_par_pre_cohfos_g1(k)]= cooccurrence_vec( Sx_sphfo_g1,Ex_sphfo_g1, Sx_pre_vec',Ex_pre_vec');
[~,~, count_cohfos_par_pre_g2(k), count_par_pre_cohfos_g2(k)] = cooccurrence_vec( Sx_sphfo_g2,Sx_sphfo_g2, Sx_pre_vec',Ex_pre_vec' ) ;
[~,~, count_cohfos_par_post_g1(k), count_par_post_cohfos_g1(k)]= cooccurrence_vec( Sx_sphfo_g1,Ex_sphfo_g1, Sx_post_vec',Ex_post_vec' );
[~,~, count_cohfos_par_post_g2(k), count_par_post_cohfos_g2(k)] = cooccurrence_vec( Sx_sphfo_g2,Sx_sphfo_g2, Sx_post_vec', Ex_post_vec' ) ;

%% Coocur shuffled PAR spindle and hfos.

parfor r=1:1000
   
    [~,~,count_cohfos_par_rand_g1(k,r), count_par_cohfos_rand_g1(k,r)]= cooccurrence_vec( Sx_sphfo_g1, Ex_sphfo_g1, Sr_par(r,:)', Er_par(r,:)' );
    [~,~, count_cohfos_par_rand_g2(k,r), count_par_cohfos_rand_g2(k,r)] = cooccurrence_vec( Sx_sphfo_g2,Ex_sphfo_g2, Sr_par(r,:)', Er_par(r,:)' );
end


%% Coocur PAR spindle and HPC ripples

[~,~,count_hpc_par(k), count_par_hpc(k)] = cooccurrence_vec(ripples_nhpc_vec(:,1),ripples_nhpc_vec(:,3), spindles_ppc_vec(:,1), spindles_ppc_vec(:,3));
[~,~,count_hpc_par_pre(k),count_par_pre_hpc(k)]=cooccurrence_vec(ripples_nhpc_vec(:,1), ripples_nhpc_vec(:,3), Sx_pre_vec',Ex_pre_vec' );
[~,~,count_hpc_par_post(k),count_par_post_hpc(k)]=cooccurrence_vec(ripples_nhpc_vec(:,1), ripples_nhpc_vec(:,3), Sx_post_vec',Ex_post_vec' );

%% Coocur PAR spindle and shuffled HPC ripples.

 parfor r=1:1000
     
     [~,~, count_hpc_par_rand(k,r), count_par_hpc_rand(k,r)] = cooccurrence_vec(ripples_nhpc_vec(:,1),ripples_nhpc_vec(:,3), Sr_par(r,:)', Er_par(r,:)' ) ;
 end

progress_bar(k,length(g),f)
    cd ..    
    end
toc

%% P values

%% NORMALIZATION (Z-SCORE-LIKE)

if(Rat==24)
    
    save('cohfos_cooccurrences_24.mat','count_cohfos_par_g1', 'count_cohfos_par_g2', ...
        'count_cohfos_par_rand_g1', 'count_cohfos_par_rand_g2', 'count_cohfos_pfc_g1', 'count_cohfos_pfc_g2', ...
        'count_cohfos_pfc_rand_g1', 'count_cohfos_pfc_rand_g2', 'count_par_cohfos_g1',  'count_par_cohfos_g2', ...
        'count_par_cohfos_rand_g1', 'count_par_cohfos_rand_g2', 'count_pfc_cohfos_g1',  'count_pfc_cohfos_g2', ...
        'count_pfc_cohfos_rand_g1', 'count_pfc_cohfos_rand_g2', 'count_hpc_par', 'count_hpc_par_rand', 'count_hpc_pfc', ...
        'count_hpc_pfc_rand', 'count_par_hpc', 'count_par_hpc_rand', 'count_pfc_hpc', 'count_pfc_hpc_rand');
    
    % PAR spindles-HPC all

    aver_hpc_par_24=count_hpc_par_rand-count_hpc_par.';
    aver_hpc_par_24=aver_hpc_par_24./std(aver_hpc_par_24.').';

    % PAR spindles-HPC Unique
    
    aver_par_hpc_24=count_par_hpc_rand-count_par_hpc.';
    aver_par_hpc_24=aver_par_hpc_24./std(aver_par_hpc_24.').';
    
    % PFC spindles-HPC all
    
    aver_hpc_pfc_24=count_hpc_pfc_rand-count_hpc_pfc.';
    aver_hpc_pfc_24=aver_hpc_pfc_24./std(aver_hpc_pfc_24.').';
    
    % PFC spindles-HPC Unique 
    
    aver_pfc_hpc_24=count_pfc_hpc_rand-count_pfc_hpc.';
    aver_pfc_hpc_24=aver_pfc_hpc_24./std(aver_pfc_hpc_24.').';
    
    % PAR spindles - slow HFOS
    
    aver_par_g1_24=count_cohfos_par_rand_g1- count_cohfos_par_g1.';
    aver_par_g1_24=aver_par_g1_24./std(aver_par_g1_24.').';
    
    % PAR spindles - fast HFOS
    
    aver_par_g2_24=count_cohfos_par_rand_g2- count_cohfos_par_g2.';
    aver_par_g2_24=aver_par_g2_24./std(aver_par_g2_24.').';
    
    % PFC spindles - slow HFOS
    
    aver_pfc_g1_24=count_cohfos_pfc_rand_g1- count_cohfos_pfc_g1.';
    aver_pfc_g1_24=aver_pfc_g1_24./std(aver_pfc_g1_24.').';
    
    % PFC spindles - fast HFOS

    aver_pfc_g2_24=count_cohfos_pfc_rand_g2- count_cohfos_pfc_g2.';
    aver_pfc_g2_24=aver_pfc_g2_24./std(aver_pfc_g2_24.').';
    
    % PFC spindle pre - slow HFOS
    
    aver_pfc_pre_g1_24=count_cohfos_pfc_rand_g1- count_cohfos_pfc_pre_g1.';
    aver_pfc_pre_g1_24=aver_pfc_pre_g1_24./std(aver_pfc_pre_g1_24.').';
    
    % PFC spindle post - slow HFOS
    
    aver_pfc_post_g1_24=count_cohfos_pfc_rand_g1- count_cohfos_pfc_post_g1.';
    aver_pfc_post_g1_24=aver_pfc_post_g1_24./std(aver_pfc_post_g1_24.').';
     
    % PFC spindle pre - fast HFOS
    
    aver_pfc_pre_g2_24=count_cohfos_pfc_rand_g2- count_cohfos_pfc_pre_g2.';
    aver_pfc_pre_g2_24=aver_pfc_pre_g2_24./std(aver_pfc_pre_g2_24.').';
    
    % PFC spindle post - fast HFOS
    
    aver_pfc_post_g2_24=count_cohfos_pfc_rand_g2- count_cohfos_pfc_post_g2.';
    aver_pfc_post_g2_24=aver_pfc_post_g2_24./std(aver_pfc_post_g2_24.').';
    
    % PAR spindle pre - slow HFOS
    
    aver_par_pre_g1_24=count_cohfos_par_rand_g1- count_cohfos_par_pre_g1.';
    aver_par_pre_g1_24=aver_par_pre_g1_24./std(aver_par_pre_g1_24.').';
    
    % PAR spindle post - slow HFOS
    
    aver_par_post_g1_24=count_cohfos_par_rand_g1- count_cohfos_par_post_g1.';
    aver_par_post_g1_24=aver_par_post_g1_24./std(aver_par_post_g1_24.').';
    
    % PAR spindle pre - fast HFOS
    
    aver_par_pre_g2_24=count_cohfos_par_rand_g2- count_cohfos_par_pre_g2.';
    aver_par_pre_g2_24=aver_par_pre_g2_24./std(aver_par_pre_g2_24.').';
        
    % PAR spindle post - fast HFOS
    
    aver_par_post_g2_24=count_cohfos_par_rand_g2- count_cohfos_par_post_g2.';
    aver_par_post_g2_24=aver_par_post_g2_24./std(aver_par_post_g2_24.').';
    
    % PFC spindle pre - ripples (unique)
    
    aver_pfc_pre_hpc_24=count_pfc_hpc_rand-count_pfc_pre_hpc.';
    aver_pfc_pre_hpc_24=aver_pfc_pre_hpc_24./std(aver_pfc_pre_hpc_24.').';
    
    % PFC spindle post - ripples
    
    aver_pfc_post_hpc_24=count_pfc_hpc_rand-count_pfc_post_hpc.';
    aver_pfc_post_hpc_24=aver_pfc_post_hpc_24./std(aver_pfc_post_hpc_24.').';
    
    % PAR spindle pre - ripples
    
    aver_par_pre_hpc_24=count_par_hpc_rand-count_par_pre_hpc.';
    aver_par_pre_hpc_24=aver_par_pre_hpc_24./std(aver_par_pre_hpc_24.').';
    
    % PAR spindle post - ripples
    
    aver_par_post_hpc_24=count_par_hpc_rand-count_par_post_hpc.';
    aver_par_post_hpc_24=aver_par_post_hpc_24./std(aver_par_post_hpc_24.').';
    
    save('aver_norm_24.mat','aver_hpc_par_24', 'aver_par_hpc_24', 'aver_hpc_pfc_24', 'aver_pfc_hpc_24', ...
        'aver_par_g1_24', 'aver_par_g2_24', 'aver_pfc_g1_24','aver_pfc_g2_24' )
    
elseif(Rat==27)
%%
    save('cohfos_cooccurrences_27.mat','count_cohfos_par_g1', 'count_cohfos_par_g2', ...
        'count_cohfos_par_rand_g1', 'count_cohfos_par_rand_g2', 'count_cohfos_pfc_g1', 'count_cohfos_pfc_g2', ...
        'count_cohfos_pfc_rand_g1', 'count_cohfos_pfc_rand_g2', 'count_par_cohfos_g1',  'count_par_cohfos_g2', ...
        'count_par_cohfos_rand_g1', 'count_par_cohfos_rand_g2', 'count_pfc_cohfos_g1',  'count_pfc_cohfos_g2', ...
        'count_pfc_cohfos_rand_g1', 'count_pfc_cohfos_rand_g2', 'count_hpc_par', 'count_hpc_par_rand', 'count_hpc_pfc', ...
        'count_hpc_pfc_rand', 'count_par_hpc', 'count_par_hpc_rand', 'count_pfc_hpc', 'count_pfc_hpc_rand' );
      
    % PAR spindles-HPC all

    aver_hpc_par_27=count_hpc_par_rand-count_hpc_par.';
    aver_hpc_par_27=aver_hpc_par_27./std(aver_hpc_par_27.').';

    % PAR spindles-HPC Unique
    
    aver_par_hpc_27=count_par_hpc_rand-count_par_hpc.';
    aver_par_hpc_27=aver_par_hpc_27./std(aver_par_hpc_27.').';
    
    % PFC spindles-HPC all
    
    aver_hpc_pfc_27=count_hpc_pfc_rand-count_hpc_pfc.';
    aver_hpc_pfc_27=aver_hpc_pfc_27./std(aver_hpc_pfc_27.').';
    
    % PFC spindles-HPC Unique 
    
    aver_pfc_hpc_27=count_pfc_hpc_rand-count_pfc_hpc.';
    aver_pfc_hpc_27=aver_pfc_hpc_27./std(aver_pfc_hpc_27.').';
    
    % PAR spindles - slow HFOS
    
    aver_par_g1_27=count_cohfos_par_rand_g1- count_cohfos_par_g1.';
    aver_par_g1_27=aver_par_g1_27./std(aver_par_g1_27.').';
    
    % PAR spindles - fast HFOS
    
    aver_par_g2_27=count_cohfos_par_rand_g2- count_cohfos_par_g2.';
    aver_par_g2_27=aver_par_g2_27./std(aver_par_g2_27.').';
    
    % PFC spindles - slow HFOS
    
    aver_pfc_g1_27=count_cohfos_pfc_rand_g1- count_cohfos_pfc_g1.';
    aver_pfc_g1_27=aver_pfc_g1_27./std(aver_pfc_g1_27.').';
    
    % PFC spindles - fast HFOS

    aver_pfc_g2_27=count_cohfos_pfc_rand_g2- count_cohfos_pfc_g2.';
    aver_pfc_g2_27=aver_pfc_g2_27./std(aver_pfc_g2_27.').';
    
    % PFC spindle pre - slow HFOS
    
    aver_pfc_pre_g1_27=count_cohfos_pfc_rand_g1- count_cohfos_pfc_pre_g1.';
    aver_pfc_pre_g1_27=aver_pfc_pre_g1_27./std(aver_pfc_pre_g1_27.').';
    
    % PFC spindle post - slow HFOS
    
    aver_pfc_post_g1_27=count_cohfos_pfc_rand_g1- count_cohfos_pfc_post_g1.';
    aver_pfc_post_g1_27=aver_pfc_post_g1_27./std(aver_pfc_post_g1_27.').';
     
    % PFC spindle pre - fast HFOS
    
    aver_pfc_pre_g2_27=count_cohfos_pfc_rand_g2- count_cohfos_pfc_pre_g2.';
    aver_pfc_pre_g2_27=aver_pfc_pre_g2_27./std(aver_pfc_pre_g2_27.').';
    
    % PFC spindle post - fast HFOS
    
    aver_pfc_post_g2_27=count_cohfos_pfc_rand_g2- count_cohfos_pfc_post_g2.';
    aver_pfc_post_g2_27=aver_pfc_post_g2_27./std(aver_pfc_post_g2_27.').';
    
    % PAR spindle pre - slow HFOS
    
    aver_par_pre_g1_27=count_cohfos_par_rand_g1- count_cohfos_par_pre_g1.';
    aver_par_pre_g1_27=aver_par_pre_g1_27./std(aver_par_pre_g1_27.').';
    
    % PAR spindle post - slow HFOS
    
    aver_par_post_g1_27=count_cohfos_par_rand_g1- count_cohfos_par_post_g1.';
    aver_par_post_g1_27=aver_par_post_g1_27./std(aver_par_post_g1_27.').';
    
    % PAR spindle pre - fast HFOS
    
    aver_par_pre_g2_27=count_cohfos_par_rand_g2- count_cohfos_par_pre_g2.';
    aver_par_pre_g2_27=aver_par_pre_g2_27./std(aver_par_pre_g2_27.').';
        
    % PAR spindle post - fast HFOS
    
    aver_par_post_g2_27=count_cohfos_par_rand_g2- count_cohfos_par_post_g2.';
    aver_par_post_g2_27=aver_par_post_g2_27./std(aver_par_post_g2_27.').';
    
    % PFC spindle pre - ripples (unique)
    
    aver_pfc_pre_hpc_27=count_pfc_hpc_rand-count_pfc_pre_hpc.';
    aver_pfc_pre_hpc_27=aver_pfc_pre_hpc_27./std(aver_pfc_pre_hpc_27.').';
    
    % PFC spindle post - ripples
    
    aver_pfc_post_hpc_27=count_pfc_hpc_rand-count_pfc_post_hpc.';
    aver_pfc_post_hpc_27=aver_pfc_post_hpc_27./std(aver_pfc_post_hpc_27.').';
    
    % PAR spindle pre - ripples
    
    aver_par_pre_hpc_27=count_par_hpc_rand-count_par_pre_hpc.';
    aver_par_pre_hpc_27=aver_par_pre_hpc_27./std(aver_par_pre_hpc_27.').';
    
    % PAR spindle post - ripples
    
    aver_par_post_hpc_27=count_par_hpc_rand-count_par_post_hpc.';
    aver_par_post_hpc_27=aver_par_post_hpc_27./std(aver_par_post_hpc_27.').';
    
    save('aver_norm_27.mat','aver_hpc_par_27', 'aver_par_hpc_27', 'aver_hpc_pfc_27', 'aver_pfc_hpc_27', ...
        'aver_par_g1_27', 'aver_par_g2_27', 'aver_pfc_g1_27','aver_pfc_g2_27' )
    
    elseif(Rat==26)
%%

   save('cohfos_cooccurrences_26.mat','count_cohfos_par_g1', 'count_cohfos_par_g2', ...
        'count_cohfos_par_rand_g1', 'count_cohfos_par_rand_g2', 'count_cohfos_pfc_g1', 'count_cohfos_pfc_g2', ...
        'count_cohfos_pfc_rand_g1', 'count_cohfos_pfc_rand_g2', 'count_par_cohfos_g1',  'count_par_cohfos_g2', ...
        'count_par_cohfos_rand_g1', 'count_par_cohfos_rand_g2', 'count_pfc_cohfos_g1',  'count_pfc_cohfos_g2', ...
        'count_pfc_cohfos_rand_g1', 'count_pfc_cohfos_rand_g2', 'count_hpc_par', 'count_hpc_par_rand', 'count_hpc_pfc', ...
        'count_hpc_pfc_rand', 'count_par_hpc', 'count_par_hpc_rand', 'count_pfc_hpc', 'count_pfc_hpc_rand' );
          
     % PAR spindles-HPC all

    aver_hpc_par_26=count_hpc_par_rand-count_hpc_par.';
    aver_hpc_par_26=aver_hpc_par_26./std(aver_hpc_par_26.').';

    % PAR spindles-HPC Unique
    
    aver_par_hpc_26=count_par_hpc_rand-count_par_hpc.';
    aver_par_hpc_26=aver_par_hpc_26./std(aver_par_hpc_26.').';
    
    % PFC spindles-HPC all
    
    aver_hpc_pfc_26=count_hpc_pfc_rand-count_hpc_pfc.';
    aver_hpc_pfc_26=aver_hpc_pfc_26./std(aver_hpc_pfc_26.').';
    
    % PFC spindles-HPC Unique 
    
    aver_pfc_hpc_26=count_pfc_hpc_rand-count_pfc_hpc.';
    aver_pfc_hpc_26=aver_pfc_hpc_26./std(aver_pfc_hpc_26.').';
    
    % PAR spindles - slow HFOS
    
    aver_par_g1_26=count_cohfos_par_rand_g1- count_cohfos_par_g1.';
    aver_par_g1_26=aver_par_g1_26./std(aver_par_g1_26.').';
    
    % PAR spindles - fast HFOS
    
    aver_par_g2_26=count_cohfos_par_rand_g2- count_cohfos_par_g2.';
    aver_par_g2_26=aver_par_g2_26./std(aver_par_g2_26.').';
    
    % PFC spindles - slow HFOS
    
    aver_pfc_g1_26=count_cohfos_pfc_rand_g1- count_cohfos_pfc_g1.';
    aver_pfc_g1_26=aver_pfc_g1_26./std(aver_pfc_g1_26.').';
    
    % PFC spindles - fast HFOS

    aver_pfc_g2_26=count_cohfos_pfc_rand_g2- count_cohfos_pfc_g2.';
    aver_pfc_g2_26=aver_pfc_g2_26./std(aver_pfc_g2_26.').';
    
    % PFC spindle pre - slow HFOS
    
    aver_pfc_pre_g1_26=count_cohfos_pfc_rand_g1- count_cohfos_pfc_pre_g1.';
    aver_pfc_pre_g1_26=aver_pfc_pre_g1_26./std(aver_pfc_pre_g1_26.').';
    
    % PFC spindle post - slow HFOS
    
    aver_pfc_post_g1_26=count_cohfos_pfc_rand_g1- count_cohfos_pfc_post_g1.';
    aver_pfc_post_g1_26=aver_pfc_post_g1_26./std(aver_pfc_post_g1_26.').';
     
    % PFC spindle pre - fast HFOS
    
    aver_pfc_pre_g2_26=count_cohfos_pfc_rand_g2- count_cohfos_pfc_pre_g2.';
    aver_pfc_pre_g2_26=aver_pfc_pre_g2_26./std(aver_pfc_pre_g2_26.').';
    
    % PFC spindle post - fast HFOS
    
    aver_pfc_post_g2_26=count_cohfos_pfc_rand_g2- count_cohfos_pfc_post_g2.';
    aver_pfc_post_g2_26=aver_pfc_post_g2_26./std(aver_pfc_post_g2_26.').';
    
    % PAR spindle pre - slow HFOS
    
    aver_par_pre_g1_26=count_cohfos_par_rand_g1- count_cohfos_par_pre_g1.';
    aver_par_pre_g1_26=aver_par_pre_g1_26./std(aver_par_pre_g1_26.').';
    
    % PAR spindle post - slow HFOS
    
    aver_par_post_g1_26=count_cohfos_par_rand_g1- count_cohfos_par_post_g1.';
    aver_par_post_g1_26=aver_par_post_g1_26./std(aver_par_post_g1_26.').';
    
    % PAR spindle pre - fast HFOS
    
    aver_par_pre_g2_26=count_cohfos_par_rand_g2- count_cohfos_par_pre_g2.';
    aver_par_pre_g2_26=aver_par_pre_g2_26./std(aver_par_pre_g2_26.').';
        
    % PAR spindle post - fast HFOS
    
    aver_par_post_g2_26=count_cohfos_par_rand_g2- count_cohfos_par_post_g2.';
    aver_par_post_g2_26=aver_par_post_g2_26./std(aver_par_post_g2_26.').';
    
    % PFC spindle pre - ripples (unique)
    
    aver_pfc_pre_hpc_26=count_pfc_hpc_rand-count_pfc_pre_hpc.';
    aver_pfc_pre_hpc_26=aver_pfc_pre_hpc_26./std(aver_pfc_pre_hpc_26.').';
    
    % PFC spindle post - ripples
    
    aver_pfc_post_hpc_26=count_pfc_hpc_rand-count_pfc_post_hpc.';
    aver_pfc_post_hpc_26=aver_pfc_post_hpc_26./std(aver_pfc_post_hpc_26.').';
    
    % PAR spindle pre - ripples
    
    aver_par_pre_hpc_26=count_par_hpc_rand-count_par_pre_hpc.';
    aver_par_pre_hpc_26=aver_par_pre_hpc_26./std(aver_par_pre_hpc_26.').';
    
    % PAR spindle post - ripples
    
    aver_par_post_hpc_26=count_par_hpc_rand-count_par_post_hpc.';
    aver_par_post_hpc_26=aver_par_post_hpc_26./std(aver_par_post_hpc_26.').';
    
    
    save('aver_norm_26.mat','aver_hpc_par_26', 'aver_par_hpc_26', 'aver_hpc_pfc_26', 'aver_pfc_hpc_26', ...
        'aver_par_g1_26', 'aver_par_g2_26', 'aver_pfc_g1_26','aver_pfc_g2_26' )
end

xo
rat = [24 26 27];

%% Tables with p-values per condition

for l=1:3
    for n=1:4
    vec = eval(['aver_hpc_par_', num2str(rat(l)), '(', num2str(n), ',:)']);
    vec2 = eval(['aver_par_hpc_', num2str(rat(l)), '(', num2str(n), ',:)']);    
    vec3 = eval(['aver_hpc_pfc_', num2str(rat(l)), '(', num2str(n), ',:)']);    
    vec4 = eval(['aver_pfc_hpc_', num2str(rat(l)), '(', num2str(n), ',:)']);        
    vec5 = eval(['aver_par_g1_', num2str(rat(l)), '(', num2str(n), ',:)']);
    vec6 = eval(['aver_par_g2_', num2str(rat(l)), '(', num2str(n), ',:)']);        
    vec7 = eval(['aver_pfc_g1_', num2str(rat(l)), '(', num2str(n), ',:)']);  
    vec8 = eval(['aver_pfc_g2_', num2str(rat(l)), '(', num2str(n), ',:)']);  
    vec9 = eval(['aver_par_pre_g1_', num2str(rat(l)), '(', num2str(n), ',:)']);
    vec10 = eval(['aver_par_pre_g2_', num2str(rat(l)), '(', num2str(n), ',:)']);        
    vec11 = eval(['aver_pfc_pre_g1_', num2str(rat(l)), '(', num2str(n), ',:)']);  
    vec12 = eval(['aver_pfc_pre_g2_', num2str(rat(l)), '(', num2str(n), ',:)']);  
    vec13 = eval(['aver_par_post_g1_', num2str(rat(l)), '(', num2str(n), ',:)']);
    vec14 = eval(['aver_par_post_g2_', num2str(rat(l)), '(', num2str(n), ',:)']);        
    vec15 = eval(['aver_pfc_post_g1_', num2str(rat(l)), '(', num2str(n), ',:)']);  
    vec16 = eval(['aver_pfc_post_g2_', num2str(rat(l)), '(', num2str(n), ',:)']);  
    vec17 = eval(['aver_par_pre_hpc_', num2str(rat(l)), '(', num2str(n), ',:)']);
    vec18 = eval(['aver_par_post_hpc_', num2str(rat(l)), '(', num2str(n), ',:)']);        
    vec19 = eval(['aver_pfc_pre_hpc_', num2str(rat(l)), '(', num2str(n), ',:)']);  
    vec20 = eval(['aver_pfc_post_hpc_', num2str(rat(l)), '(', num2str(n), ',:)']);  
    
    
    
       
    pv_individual(l,n)=(1+sum(vec >=0))/(length(vec)+1);
    pv2_individual(l,n)=(1+sum(vec2 >=0))/(length(vec2)+1);
    pv3_individual(l,n)=(1+sum(vec3 >=0))/(length(vec3)+1);
    pv4_individual(l,n)=(1+sum(vec4 >=0))/(length(vec4)+1);
    pv5_individual(l,n)=(1+sum(vec5 >=0))/(length(vec5)+1);
    pv6_individual(l,n)=(1+sum(vec6 >=0))/(length(vec6)+1);
    pv7_individual(l,n)=(1+sum(vec7 >=0))/(length(vec7)+1);
    pv8_individual(l,n)=(1+sum(vec8 >=0))/(length(vec8)+1);
    pv9_individual(l,n)=(1+sum(vec9 >=0))/(length(vec9)+1);
    pv10_individual(l,n)=(1+sum(vec10 >=0))/(length(vec10)+1);
    pv11_individual(l,n)=(1+sum(vec11 >=0))/(length(vec11)+1);
    pv12_individual(l,n)=(1+sum(vec12 >=0))/(length(vec12)+1);
    pv13_individual(l,n)=(1+sum(vec13 >=0))/(length(vec13)+1);
    pv14_individual(l,n)=(1+sum(vec14 >=0))/(length(vec14)+1);
    pv15_individual(l,n)=(1+sum(vec15 >=0))/(length(vec15)+1);
    pv16_individual(l,n)=(1+sum(vec16 >=0))/(length(vec16)+1);
    pv17_individual(l,n)=(1+sum(vec17 >=0))/(length(vec17)+1);
    pv18_individual(l,n)=(1+sum(vec18 >=0))/(length(vec18)+1);
    pv19_individual(l,n)=(1+sum(vec19 >=0))/(length(vec19)+1);
    pv20_individual(l,n)=(1+sum(vec20 >=0))/(length(vec20)+1);
    

    end
end

%% TABLE
cd ..

%P values of individual conditions and animal
% 1. PAR - HPC all
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A1')
            writematrix('Coocur PAR HPC all', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A32')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',2,'Range','B33:L37')   

% 2. PAR - HPC unique
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv2_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A1')
            writematrix('Coocur PAR HPC unique', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A37')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',2,'Range','B38:L42')   
            
% 3. PFC - HPC all
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv3_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g ,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A1')
            writematrix('Coocur PFC - HPC all', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A2')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',2,'Range','B3:L7')   
            
% 4. PFC - HPC unique
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv4_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A1')
            writematrix('Coocur PFC - HPC unique', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A7')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',2,'Range','B8:L12')   
            
% 5. PAR - slow HFOS
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv5_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 1, 'Range', 'A1')
            writematrix('Coocur PAR - slow HFOS ', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A12')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',2,'Range','B13:L17')   

% 6. PAR - fast HFOS
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv6_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A1')
            writematrix('Coocur PAR - fast HFOs', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A17')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',2,'Range','B18:L22')   

% 7. PFC - slow HFOS
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv7_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A1')
            writematrix('Coocur PFC - slow HFOS', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A22')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',2,'Range','B23:L27')   

% 8. PFC - fast HFOS
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv8_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A1')
            writematrix('Coocur PFC - fast HFOS', 'P_values_Lisa_70.xls','Sheet', 2, 'Range', 'A27')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',2,'Range','B28:L32')   
            
%% pre spindles

% 9. PFC pre - slow HFOS
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv11_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A1')
            writematrix('Coocur PFC pre - Slow HFOS', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A2')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',3,'Range','B3:L7')   

% 10. PFC pre - fast HFOS
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv12_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A1')
            writematrix('Coocur PFC pre - fast HFOS', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A7')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',3,'Range','B8:L12')   

% 13. PAR pre - slow HFOS
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv9_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A1')
            writematrix('Coocur PAR pre - slow HFOS', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A22')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',3,'Range','B23:L27')   

% 14. PAR pre - fast HFOS
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv10_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A1')
            writematrix('Coocur PAR pre - fast HFOS', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A27')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',3,'Range','B28:L32')   

% 17. PFC -pre HPC unique
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv19_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A1')
            writematrix('Coocur PFC pre - HPC unique', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A44')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',3,'Range','B43:L47')   

% 19. PAR pre - HPC unique
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv17_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A1')
            writematrix('Coocur PAR pre - HPC unique', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A52')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',3,'Range','B53:L57')  
           
%% post spindles

% 11. PFC post - slow HFOS
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv15_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A1')
            writematrix('Coocur PFC post - slow HFOS', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A12')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',3,'Range','B13:L17')   
            
% 12. PFC post - fast HFOS
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv16_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A1')
            writematrix('Coocur PFC post - fast HFOS', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A17')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',3,'Range','B18:L22')   

% 15. PAR post - slow HFOS
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv13_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A1')
            writematrix('Coocur PAR post - slow HFOS', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A32')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',3,'Range','B33:L37')   
            
% 16. PAR post - fast HFOS
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv14_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A1')
            writematrix('Coocur PAR post - fast HFOS', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A37')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',3,'Range','B38:L42')   
            
% 18. PFC post- HPC unique
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv20_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A1')
            writematrix('Coocur PFC -post HPC unique', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A47')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',3,'Range','B48:L52')               
                        
% 20. PAR post- HPC unique
    TT=table;
    TT.Variables= [[{'Rat 24'}; {'Rat 26'}; {'Rat 27'}] num2cell(pv18_individual)];
    
    TT.Properties.VariableNames=['Rat' cellfun(@(equis) equis,g,'UniformOutput',false)];
            
            writematrix('P values all conditions all animals', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A1')
            writematrix('Coocur PAR post - HPC unique', 'P_values_Lisa_70.xls','Sheet', 3, 'Range', 'A57')
            writetable(TT,'P_values_Lisa_70.xls','Sheet',3,'Range','B58:L62')   
 