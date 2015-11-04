function step3_cluster(opts,data)

timerA=tic;

fprintf('Step 3: Cluster... ');

detect_times_prefix=opts.detect_times_prefix;
detect_clips_prefix=opts.detect_clips_prefix;
cluster_labels_prefix=opts.cluster_labels_prefix;
cluster_waveforms_prefix=opts.cluster_waveforms_prefix;
cluster_features_prefix=opts.cluster_features_prefix;

X=data.X;

AM=readmda(opts.adjacency);

fprintf('Clustering...\n');
for j=1:size(X,1)
    fname_detect_times=[detect_times_prefix,sprintf('%d.mda',j)];
    fname_detect_clips=[detect_clips_prefix,sprintf('%d.mda',j)];
    fname_cluster_labels=[cluster_labels_prefix,sprintf('%d.mda',j)];
    fname_cluster_waveforms=[cluster_waveforms_prefix,sprintf('%d.mda',j)];
    fname_cluster_features=[cluster_features_prefix,sprintf('%d.mda',j)];
    if (~exist(fname_cluster_labels,'file'))||(~exist(fname_cluster_waveforms,'file'))||(~exist(fname_cluster_features,'file'))
        fprintf('Reading %s... ',fname_detect_times);
        detect_times=readmda(fname_detect_times);
        fprintf('Reading %s... ',fname_detect_clips);
        clips=readmda(fname_detect_clips);
        [M,T,Nclips]=size(clips);
        fprintf('Extracting %d features...',opts.num_cluster_features);
        AM(j,j)=0;
        adjacent_channels=[j;find(AM(:,j))];
        FF=ss_eventfeatures(clips(adjacent_channels,:,:));
        FF=FF(1:opts.num_cluster_features,:);
        
        fprintf('ISO-SPLIT...\n');
        labels=isosplit(FF);
        
        fprintf('Computing average waveforms... ');
        K=max(labels);
        WF=zeros(M,T,K);
        for k=1:K
            inds=find(labels==k);
            %WF(:,:,k)=mean(clips(:,:,inds),3);
            WF(:,:,k)=median(clips(:,:,inds),3);
        end;
        fprintf('Writing %s... ',fname_cluster_waveforms);
        writemda(WF,fname_cluster_waveforms);
        fprintf('Writing %s (%d clusters detected)...\n',fname_cluster_labels,K);
        writemda(labels,fname_cluster_labels);
        fprintf('Writing %s... ',fname_cluster_features);
        writemda(FF,fname_cluster_features);
    else
        fprintf('Files exist %s, %s and %s...\n',fname_cluster_labels,fname_cluster_waveforms,fname_cluster_features);
    end;
end;

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');

end
