% plotting figures for time frequency analysis
function []=plot_granger_tf(g_tf_sp, g_tf_fp, g_tf_p, event1, event2, event3)

pos=[ 2 3; 
      1 3;
      1 2];

% frequency range is 0-300Hz 
titles = [{'PAR->PFC'} {'PAR->HPC'};{'PFC->PAR'} {'PFC->HPC'}; {'HPC->PAR'} {'HPC->PFC'}];
for i=1:3
    for j=1:2
        figure,
        %tf_p=squeeze(g_tf_p.grangerspctrm(i,pos(i,j),:,:));
        tf_p=mean(squeeze(g_tf_p.grangerspctrm(i,pos(i,j),:,:,:)), 3);
        tf_p=tf_p';
        
        %tf_sp=squeeze(g_tf_sp.grangerspctrm(i,pos(i,j),:,:));
        tf_sp=mean(squeeze(g_tf_sp.grangerspctrm(i,pos(i,j),:,:,:)), 3);
        tf_sp=tf_sp';
        
        %tf_fp=squeeze(g_tf_fp.grangerspctrm(i,pos(i,j),:,:));
        tf_fp=mean(squeeze(g_tf_fp.grangerspctrm(i,pos(i,j),:,:,:)), 3);
        tf_fp=tf_fp';
        
        zmin= min([min(tf_p, [],'all'), min(tf_sp, [],'all'), min(tf_fp, [],'all')],[],'all');
        zmax= max([max(tf_p, [], 'all'), max(tf_sp, [],'all'), max(tf_fp, [],'all')],[],'all');
        clim =[zmin zmax]
        subplot(1,3,1)
        imagesc(g_tf_p.time, g_tf_p.freq, tf_p, clim );
        axis xy % flip vertically
        colorbar
        colormap(hot(256))
        xlabel('time')
        ylabel('frequency')
        title([event3 '-' cell2mat(titles(i,j))])
        
        subplot(1,3,2)
        imagesc(g_tf_sp.time, g_tf_sp.freq, tf_sp, clim);
        axis xy % flip vertically
        colorbar
        colormap(hot(256))
        xlabel('time')
        ylabel('frequency')
        title([event1 '-' cell2mat(titles(i,j))])
         
        subplot(1,3,3)
        imagesc(g_tf_fp.time, g_tf_fp.freq, tf_fp, clim);
        axis xy % flip vertically
        colorbar
        colormap(hot(256))
        xlabel('time')
        ylabel('frequency')
        title([event2 '-' cell2mat(titles(i,j))])
        
    end
end

% frequency range is 0-20Hz
for i=1:3
    for j=1:2
        figure,
        tf_p=mean(squeeze(g_tf_p.grangerspctrm(i,pos(i,j),:,:,:)),3);
        tf_p=tf_p';
        tf_p20=tf_p(:,1:20);
        
        tf_sp=mean(squeeze(g_tf_sp.grangerspctrm(i,pos(i,j),:,:,:)), 3);
        tf_sp=tf_sp';
        tf_sp20=tf_sp(:,1:20);
        
        tf_fp=mean(squeeze(g_tf_fp.grangerspctrm(i,pos(i,j),:,:,:)), 3);
        tf_fp=tf_fp';
        tf_fp20=tf_fp(:,1:20);
        
        zmin= min([min(tf_p20, [],'all'), min(tf_sp20, [],'all'), min(tf_fp20, [],'all')],[],'all');
        zmax= max([max(tf_p20, [], 'all'), max(tf_sp20, [],'all'), max(tf_fp20, [],'all')],[],'all');
        clim =[zmin zmax]
        subplot(1,3,1)
        imagesc(g_tf_p.time, 1:20, tf_p20, clim );
        axis xy % flip vertically
        colorbar
        colormap(hot(256))
        xlabel('time')
        ylabel('frequency')
        title([event3 '-' cell2mat(titles(i,j))])
        
        subplot(1,3,2)
        imagesc(g_tf_sp.time, 1:20, tf_sp20, clim);
        axis xy % flip vertically
        colorbar
        colormap(hot(256))
        xlabel('time')
        ylabel('frequency')
        title([event1 '-' cell2mat(titles(i,j))])
         
        subplot(1,3,3)
        imagesc(g_tf_fp.time, 1:20, tf_fp20, clim);
        axis xy % flip vertically
        colorbar
        colormap(hot(256))
        xlabel('time')
        ylabel('frequency')
        title([event2 '-' cell2mat(titles(i,j))])
        
    end
end

freqrange=[1 300]
% permutation test
for i =1:3
    for j=1:2
        a=[i, pos(i,j)];
        zmap1=stats_high_gran(g_tf_p,g_tf_sp,freqrange, a);
        zmap2=stats_high_gran(g_tf_p,g_tf_fp,freqrange, a);
        zmap3=stats_high_gran(g_tf_sp,g_tf_fp,freqrange, a);
        
        zmin= min([min(zmap1, [],'all'), min(zmap2, [],'all'), min(zmap3, [],'all')],[],'all');
        zmax= max([max(zmap1, [], 'all'), max(zmap2, [],'all'), max(zmap3, [],'all')],[],'all');
        
        clim =[zmin zmax]
        figure,
        subplot(1, 3, 1)
        imagesc(g_tf_sp.time, g_tf_sp.freq, zmap1', clim)
        colorbar
        colormap(bluewhitered(256))
        xlabel('time')
        ylabel('frequency')
        title([event1 '(' event3 ')' '-' cell2mat(titles(i,j))])
        
        
        subplot(1, 3, 2)
        imagesc(g_tf_sp.time, g_tf_sp.freq, zmap2', clim)
        colorbar
        colormap(bluewhitered(256))
        xlabel('time')
        ylabel('frequency')
        title([event2 '(' event3 ')' '-' cell2mat(titles(i,j))])
        
        
        subplot(1, 3, 3)
        imagesc(g_tf_sp.time, g_tf_sp.freq, zmap3', clim)
        colorbar
        colormap(bluewhitered(256))
        xlabel('time')
        ylabel('frequency')
        title([event2 '(' event1 ')' '-' cell2mat(titles(i,j))])
    end
end

% permutation test
freqrange=[1 20];
for i =1:3
    for j=1:2
        a=[i, pos(i,j)];
        zmap1=stats_high_gran(g_tf_p,g_tf_sp, freqrange, a);
        zmap2=stats_high_gran(g_tf_p,g_tf_fp, freqrange, a);
        zmap3=stats_high_gran(g_tf_sp,g_tf_fp, freqrange, a);
        
        zmin= min([min(zmap1, [],'all'), min(zmap2, [],'all'), min(zmap3, [],'all')],[],'all');
        zmax= max([max(zmap1, [], 'all'), max(zmap2, [],'all'), max(zmap3, [],'all')],[],'all');
        
        clim =[zmin zmax]
        figure,
        subplot(1, 3, 1)
        imagesc(g_tf_sp.time, g_tf_sp.freq(1:20), zmap1', clim)
        colorbar
        colormap(bluewhitered(256))
        xlabel('time')
        ylabel('frequency')
        title([event1 '(' event3 ')' '-' cell2mat(titles(i,j))])
        
        
        subplot(1, 3, 2)
        imagesc(g_tf_sp.time, g_tf_sp.freq(1:20), zmap2', clim)
        colorbar
        colormap(bluewhitered(256))
        xlabel('time')
        ylabel('frequency')
        title([event2 '(' event3 ')' '-' cell2mat(titles(i,j))])
        
        
        subplot(1, 3, 3)
        imagesc(g_tf_sp.time, g_tf_sp.freq(1:20), zmap3', clim)
        colorbar
        colormap(bluewhitered(256))
        xlabel('time')
        ylabel('frequency')
        title([event2 '(' event1 ')' '-' cell2mat(titles(i,j))])
    end
end

end