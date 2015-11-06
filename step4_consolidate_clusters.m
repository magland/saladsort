function step5_consolidate_clusters(opts)

timerA=tic;

fprintf('Step 5: Consolidate clusters... ');

cluster_times_prefix=opts.cluster_times_prefix;
cluster_labels_prefix=opts.cluster_labels_prefix;
cluster_waveforms_prefix=opts.cluster_waveforms_prefix;
%cluster2_times_prefix=opts.cluster2_times_prefix;
%cluster2_labels_prefix=opts.cluster2_labels_prefix;
%cluster2_waveforms_prefix=opts.cluster2_waveforms_prefix;
consolidate_times_path=opts.consolidate_times_path;
consolidate_labels_path=opts.consolidate_labels_path;
consolidate_waveforms_path=opts.consolidate_waveforms_path;
consolidate_load_channels_path=opts.consolidate_load_channels_path;

AM=readmda(opts.adjacency);

M=size(AM,1);

TIMES=[];
LABELS=[];
WAVEFORMS=zeros(M,opts.clip_size,0);
LOAD_CHANNELS=[];

current_label=1;

for j=1:M
    fname_cluster2_times=[cluster_times_prefix,sprintf('%d.mda',j)];
    fname_cluster2_labels=[cluster_labels_prefix,sprintf('%d.mda',j)];
    fname_cluster2_waveforms=[cluster_waveforms_prefix,sprintf('%d.mda',j)];
    times=readmda(fname_cluster2_times);
    labels=readmda(fname_cluster2_labels);
    WF=readmda(fname_cluster2_waveforms);
    
    if (size(WF(:)>1))
        sizes=squeeze(sum(WF.^2,2));
        max_sizes=(max(sizes,[],1));
        rel_sizes=sizes(j,:)./max_sizes;
        labels_to_use=find(rdel_sizes>=0.9);

        for ii=1:length(labels_to_use)
            k=labels_to_use(ii);
            time_indices=find(labels==k);
            TIMES=[TIMES,times(time_indices)];
            LABELS=[LABELS,ones(size(labels(time_indices)))*current_label];
            WAVEFORMS=cat(3,WAVEFORMS,WF(:,:,k));
            LOAD_CHANNELS=[LOAD_CHANNELS,j];
            current_label=current_label+1;
        end;
    end;
end;

[TIMES,sort_inds]=sort(TIMES);
LABELS=LABELS(sort_inds);

fprintf('Writing %s...\n',consolidate_times_path);
writemda(TIMES,consolidate_times_path);
fprintf('Writing %s...\n',consolidate_labels_path);
writemda(LABELS,consolidate_labels_path);
fprintf('Writing %s...\n',consolidate_waveforms_path);
writemda(WAVEFORMS,consolidate_waveforms_path);
fprintf('Writing %s...\n',consolidate_load_channels_path);
writemda(LOAD_CHANNELS,consolidate_load_channels_path);

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');

end
