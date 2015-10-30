function step3_cluster(opts,data)

timerA=tic;

fprintf('Step 3: Cluster... ');

detect_clips_prefix=opts.detect_clips_prefix
detect_labels_prefix=opts.detect_labels_prefix;
detect_waveforms_prefix=opts.detect_waveforms_prefix;
detect_features_prefix=opts.detect_features_prefix;

X=data.X;

AM=readmda(opts.adjacency);

fprintf('Clustering...\n');
for j=1:size(X,1)
    fname_clips=[detect_clips_prefix,sprintf('%d.mda',j)];
    fname_labels=[detect_labels_prefix,sprintf('%d.mda',j)];
    fname_waveforms=[detect_waveforms_prefix,sprintf('%d.mda',j)];
    fname_features=[detect_features_prefix,sprintf('%d.mda',j)];
    if ((~exist(fname_labels,'file'))||(~exist(fname_waveforms,'file'))||(~exist(fname_features,'file')))
        fprintf('Reading %s... ',fname_clips);
        clips=readmda(fname_clips);
        [M,T,Nclips]=size(clips);
        fprintf('Extracting %d features...',opts.num_cluster_features);
        adjacent_channels=find(AM(:,j));
        FF=ss_eventfeatures(clips(adjacent_channels,:,:));
        FF=FF(1:opts.num_cluster_features,:);
        fprintf('Writing %s...',fname_features);
        writemda(FF,fname_features);
        fprintf('ISO-SPLIT...');
        labels=isosplit(FF);
        fprintf('Writing %s (%d clusters detected)...\n',fname_labels,max(labels));
        writemda(labels,fname_labels);
        fprintf('Computing average waveforms... ');
        K=max(labels);
        WF=zeros(M,T,K);
        for k=1:K
            inds=find(labels==k);
            %WF(:,:,k)=mean(clips(:,:,inds),3);
            WF(:,:,k)=median(clips(:,:,inds),3);
        end;
        fprintf('Writing %s...\n',fname_waveforms);
        writemda(WF,fname_waveforms);
    else
        fprintf('Files exist %s, %s and %s...\n',fname_labels,fname_waveforms,fname_features);
    end;
end;

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');

end
