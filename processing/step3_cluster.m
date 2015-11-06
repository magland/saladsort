function step3_cluster(opts,data)

timerA=tic;

fprintf('Step 3: Cluster... ');

detect_times_prefix=opts.detect_times_prefix;
detect_clips_prefix=opts.detect_clips_prefix;
cluster_times_prefix=opts.cluster_times_prefix;
cluster_labels_prefix=opts.cluster_labels_prefix;
cluster_waveforms_prefix=opts.cluster_waveforms_prefix;
cluster_features_prefix=opts.cluster_features_prefix;

X=data.X;

AM=readmda(opts.adjacency);

fprintf('Clustering...\n');
for j=1:size(X,1)
    fname_detect_times_pos=[detect_times_prefix,sprintf('pos_%d.mda',j)];
    fname_detect_times_neg=[detect_times_prefix,sprintf('neg_%d.mda',j)];
    fname_detect_clips_pos=[detect_clips_prefix,sprintf('pos_%d.mda',j)];
    fname_detect_clips_neg=[detect_clips_prefix,sprintf('neg_%d.mda',j)];
    fname_cluster_times=[cluster_times_prefix,sprintf('%d.mda',j)];
    fname_cluster_labels=[cluster_labels_prefix,sprintf('%d.mda',j)];
    fname_cluster_labels_pos=[cluster_labels_prefix,sprintf('pos_%d.mda',j)];
    fname_cluster_labels_neg=[cluster_labels_prefix,sprintf('neg_%d.mda',j)];
    fname_cluster_waveforms=[cluster_waveforms_prefix,sprintf('%d.mda',j)];
    fname_cluster_features_pos=[cluster_features_prefix,sprintf('pos_%d.mda',j)];
    fname_cluster_features_neg=[cluster_features_prefix,sprintf('neg_%d.mda',j)];
    if (~exist(fname_cluster_labels,'file'))||(~exist(fname_cluster_waveforms,'file'))||(~exist(fname_cluster_features_pos,'file'))||(~exist(fname_cluster_features_neg,'file'))
        fprintf('Reading %s and %s...\n',fname_detect_times_pos,fname_detect_times_neg);
        detect_times_pos=readmda(fname_detect_times_pos);
        detect_times_neg=readmda(fname_detect_times_neg);
        fprintf('Reading %s and %s...\n',fname_detect_clips_pos,fname_detect_clips_neg);
        clips_pos=readmda(fname_detect_clips_pos);
        clips_neg=readmda(fname_detect_clips_neg);
        [M,T,Nclips_pos]=size(clips_pos);
        [M,T,Nclips_neg]=size(clips_neg);
        fprintf('Extracting %d features...',opts.num_cluster_features);
        AM(j,j)=0;
        adjacent_channels=[j;find(AM(:,j))];
        %mean_clips_pos=mean(clips_pos(adjacent_channels,:,:),3);
        %mean_clips_neg=mean(clips_neg(adjacent_channels,:,:),3);
        %FF_pos=ss_eventfeatures(clips_pos(adjacent_channels,:,:)-repmat(mean_clips_pos,1,1,Nclips_pos));
        %FF_neg=ss_eventfeatures(clips_neg(adjacent_channels,:,:)-repmat(mean_clips_neg,1,1,Nclips_neg));
        
        %norms_pos=repmat(sqrt(sum(sum(clips_pos(adjacent_channels,:,:).^2,1),2)),length(adjacent_channels),T,1);
        %norms_neg=repmat(sqrt(sum(sum(clips_neg(adjacent_channels,:,:).^2,1),2)),length(adjacent_channels),T,1);
        %clipsFF_pos=clips_pos(adjacent_channels,:,:)./norms_pos;
        %clipsFF_neg=clips_neg(adjacent_channels,:,:)./norms_neg;
        
        clipsFF_pos=clips_pos(adjacent_channels,:,:);
        clipsFF_neg=clips_neg(adjacent_channels,:,:);
        
        FF_pos=ss_eventfeatures(clipsFF_pos);
        FF_neg=ss_eventfeatures(clipsFF_neg);
        FF_pos=FF_pos(1:opts.num_cluster_features,:);
        FF_neg=FF_neg(1:opts.num_cluster_features,:);
        
        fprintf('ISO-SPLIT...\n');
        times=[];
        labels=[];
        if Nclips_pos>0
            labels_pos=isosplit(FF_pos);
            Kpos=max(labels_pos);
            times=[times,detect_times_pos];
            labels=[labels,labels_pos];
        else
            labels_pos=[];
            Kpos=0;
        end;
        if Nclips_neg>0
            labels_neg=isosplit(FF_neg);
            times=[times,detect_times_neg];
            labels=[labels,labels_neg+Kpos];
        else
            labels_neg=[];
        end;
        
        clips=cat(3,clips_pos,clips_neg);
        [times,sort_inds]=sort(times);
        labels=labels(sort_inds);
        clips=clips(:,:,sort_inds);
        
        fprintf('Computing average waveforms... ');
        K=max(labels);
        WF=zeros(M,T,K);
        for k=1:K
            inds=find(labels==k);
            WF(:,:,k)=median(clips(:,:,inds),3);
        end;
        
%         fprintf('Removing the noise clusters...\n');
%         noise_clusters=[];
%         for k=1:K
%             if (is_noise_cluster(
%         end;
        
%         fprintf('Removing the zero clusters...\n');
%         magnitudes=squeeze(sqrt(sum(sum(WF.^2,1),2)));
%         [~,mag_sort_inds]=sort(magnitudes);
%         zero_label_1=mag_sort_inds(1);
%         zero_label_2=mag_sort_inds(2);
%         labels(find(labels==zero_label_1))=0;
%         labels(find(labels==zero_label_2))=0;
%         okay_inds=find(labels>0);
%         times=times(okay_inds);
%         labels=labels(okay_inds);
%         clips=clips(:,:,okay_inds);
%         labels0=labels;
%         labels(find(labels0>zero_label_1))=labels(find(labels0>zero_label_1))-1;
%         labels(find(labels0>zero_label_2))=labels(find(labels0>zero_label_2))-1;
%         if (K-2~=max(labels)) error('Unexpected problem!'); end;
        
        [times,labels,clips]=remove_small_clusters(times,labels,clips,60);
        K=max(labels);
        
        WF=zeros(M,T,K);
        for k=1:K
            inds=find(labels==k);
            WF(:,:,k)=median(clips(:,:,inds),3);
        end;
        
        fprintf('Writing %s... ',fname_cluster_waveforms);
        writemda(WF,fname_cluster_waveforms);
        fprintf('Writing %s...\n',fname_cluster_times);
        writemda(times,fname_cluster_times);
        fprintf('Writing %s (%d clusters detected)...\n',fname_cluster_labels,K);
        writemda(labels,fname_cluster_labels);
        fprintf('Writing %s and %s...\n',fname_cluster_labels_pos,fname_cluster_labels_neg);
        writemda(labels_pos,fname_cluster_labels_pos);
        writemda(labels_neg,fname_cluster_labels_neg);
        fprintf('Writing %s and %s... ',fname_cluster_features_pos,fname_cluster_features_neg);
        writemda(FF_pos,fname_cluster_features_pos);
        writemda(FF_neg,fname_cluster_features_neg);
    else
        fprintf('Files exist %s, %s, %s and %s...\n',fname_cluster_labels,fname_cluster_waveforms,fname_cluster_features_pos,fname_cluster_features_neg);
    end;
end;

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');

end

function [times2,labels2,clips2]=remove_small_clusters(times,labels,clips,minsize)

K=max(labels);
mapping=zeros(1,K);
k2=0;
for k=1:K
    ct=length(find(labels==k));
    if (ct>=minsize)
        mapping(k)=k2+1; k2=k2+1;
    else
        mapping(k)=0;
    end;
end;
labels2=mapping(labels);
inds=find(labels2>0);
labels2=labels2(inds);
times2=times(inds);
clips2=clips(:,:,inds);

end

