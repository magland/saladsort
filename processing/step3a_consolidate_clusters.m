function step3a_consolidate_clusters(opts,data)

timerA=tic;

fprintf('Step 3a: Consolidate clusters... ');

detect_times_prefix=opts.detect_times_prefix;
cluster_labels_prefix=opts.cluster_labels_prefix;
cluster_waveforms_prefix=opts.cluster_waveforms_prefix;
cluster_times_path=opts.cluster_times_path;
cluster_labels_path=opts.cluster_labels_path;
cluster_waveforms_path=opts.cluster_waveforms_path;
cluster_load_channels_path=opts.cluster_load_channels_path;

AM=readmda(opts.adjacency);

X=data.X;
M=size(X,1);
N=size(X,2);

TIMES=[];
LABELS=[];
WAVEFORMS=zeros(M,opts.clip_size,0);
LOAD_CHANNELS=[];

current_label=1;

for j=1:size(X,1)
    fname_detect_times=[detect_times_prefix,sprintf('%d.mda',j)];
    fname_cluster_labels=[cluster_labels_prefix,sprintf('%d.mda',j)];
    fname_cluster_waveforms=[cluster_waveforms_prefix,sprintf('%d.mda',j)];
    times=readmda(fname_detect_times);
    labels=readmda(fname_cluster_labels);
    WF=readmda(fname_cluster_waveforms);
    
    sizes=squeeze(sum(WF.^2,2));
    max_sizes=(max(sizes,[],1));
    rel_sizes=sizes(j,:)./max_sizes;
    labels_to_use=find(rel_sizes>=0.9);
    
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

[TIMES,sort_inds]=sort(TIMES);
LABELS=LABELS(sort_inds);

fprintf('Writing %s...\n',cluster_times_path);
writemda(TIMES,cluster_times_path);
fprintf('Writing %s...\n',cluster_labels_path);
writemda(LABELS,cluster_labels_path);
fprintf('Writing %s...\n',cluster_waveforms_path);
writemda(WAVEFORMS,cluster_waveforms_path);
fprintf('Writing %s...\n',cluster_load_channels_path);
writemda(LOAD_CHANNELS,cluster_load_channels_path);

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');

end
