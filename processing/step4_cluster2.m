function step4_cluster2(opts,data)

timerA=tic;

fprintf('Step 4: Cluster2... ');

detect_times_prefix=opts.detect_times_prefix;
cluster_waveforms_prefix=opts.cluster_waveforms_prefix;
cluster2_sliding_ips_prefix=opts.cluster2_sliding_ips_prefix;
cluster2_times_prefix=opts.cluster2_times_prefix;
cluster2_labels_prefix=opts.cluster2_labels_prefix;
cluster2_waveforms_prefix=opts.cluster2_waveforms_prefix;
cluster2_features_prefix=opts.cluster2_features_prefix;

X=data.X;
N=size(X,2);

AM=readmda(opts.adjacency);

fprintf('Clustering...\n');
for j=1:size(X,1)
    fname_detect_times=[detect_times_prefix,sprintf('%d.mda',j)];
    fname_cluster_waveforms=[cluster_waveforms_prefix,sprintf('%d.mda',j)];
    fname_cluster2_sliding_ips=[cluster2_sliding_ips_prefix,sprintf('%d.mda',j)];
    fname_cluster2_times=[cluster2_times_prefix,sprintf('%d.mda',j)];
    fname_cluster2_labels=[cluster2_labels_prefix,sprintf('%d.mda',j)];
    fname_cluster2_waveforms=[cluster2_waveforms_prefix,sprintf('%d.mda',j)];
    fname_cluster2_features=[cluster2_features_prefix,sprintf('%d.mda',j)];
    if ((~exist(fname_cluster2_sliding_ips,'file'))||(~exist(fname_cluster2_times,'file'))||(~exist(fname_cluster2_labels,'file'))||(~exist(fname_cluster2_waveforms,'file'))||(~exist(fname_cluster2_features,'file')))
        fprintf('Reading %s... ',fname_detect_times);
        detect_times=readmda(fname_detect_times);
        fprintf('Reading %s... ',fname_cluster_waveforms);
        cluster_waveforms=readmda(fname_cluster_waveforms);
        
        fprintf('Identifying clusters to use...\n');
        sizes=squeeze(sum(cluster_waveforms.^2,2));
        max_sizes=(max(sizes,[],1));
        rel_sizes=sizes(j,:)./max_sizes;
        labels_to_use=find(rel_sizes>=0.9);
        
        fprintf('Auto truncating waveforms...\n');
        WF0=cluster_waveforms(j,:,labels_to_use); %only the main channel!
        K1=size(WF0,3);
        maxval=max(abs(WF0(:)));
        tmp=sum(abs(WF0)>=maxval*0.1,3);
        inds00=find(tmp~=0);
        inds00min=inds00(1); inds00max=inds00(end);
        WF=WF0(1,inds00min:inds00max,:);
        fprintf('%d x %d x %d\n',size(WF,1),size(WF,2),size(WF,3));
        T=size(WF,2);
        
        fprintf('Computing sliding ips... ');
        sliding_ips=zeros(K1,N);
        tt0=-ceil((T+1)/2); tt1=tt0+T-1;
        disp([tt0,tt1]);
        aa=(-tt0+1):(N-tt1);
        for k=1:K1
            fprintf('|');
            tmp=zeros(1,N);
            for dt=tt0:tt1
                fprintf('.');
                tmp(aa)=tmp(aa)+WF0(1,dt-tt0+1,k)*X(j,aa+dt);
            end;
            tmp=tmp/sqrt(var(tmp));
            sliding_ips(k,:)=tmp;
        end;
        fprintf('Writing %s...\n',fname_cluster2_sliding_ips);
        writemda(sliding_ips,fname_cluster2_sliding_ips);
        
        fprintf('Detecting events... ');
        sliding_ips_maxabs=max(abs(sliding_ips),[],1);
        critical_times=find(abs(sliding_ips_maxabs)>=4);
        %Now make sure we only use those that are global maxima over the radius
        %T/2
        T_2=ceil((T+1)/2);
        use_it=ones(1,length(critical_times));
        for ii=1:length(critical_times)
            if ((critical_times(ii)-T_2)<1) use_it(ii)=0; end;
            if ((critical_times(ii)+T_2)>N) use_it(ii)=0; end;
            k=ii-1;
            while (k>=1)&&(critical_times(k)>=critical_times(ii)-T_2)
                if (sliding_ips_maxabs(critical_times(k))>sliding_ips_maxabs(critical_times(ii)))
                    use_it(ii)=0;
                else
                    use_it(k)=0;
                end;
                k=k-1;
            end
        end;
        critical_times=critical_times(find(use_it));
        if (length(critical_times)==0) error('Unexpected (but not impossible) problem... found no critical times'); end;
        fprintf('Writing %s...\n',fname_cluster2_times);
        writemda(critical_times,fname_cluster2_times);
        
        fprintf('Extracting clips... ');
        adjacent_channels=find(AM(:,j));
        clips=extract_clips(X(adjacent_channels,:),critical_times,T);
        [M,T,Nclips]=size(clips);
        fprintf('Extracting %d features...',opts.num_cluster_features);
        FF=ss_eventfeatures(clips);
        FF=FF(1:opts.num_cluster_features,:);
        fprintf('Writing %s...\n',fname_cluster2_features);
        writemda(FF,fname_cluster2_features);
        
        fprintf('ISO-SPLIT...\n');
        labels=isosplit(FF);
        fprintf('Writing %s (%d clusters detected)...\n',fname_cluster2_labels,max(labels));
        writemda(labels,fname_cluster2_labels);
        
        fprintf('Computing average waveforms... ');
        K=max(labels);
        WF0=zeros(M,T,K);
        for k=1:K
            inds=find(labels==k);
            %WF(:,:,k)=mean(clips(:,:,inds),3);
            WF0(:,:,k)=median(clips(:,:,inds),3);
        end;
        fprintf('Writing %s... ',fname_cluster2_waveforms);
        writemda(WF0,fname_cluster2_waveforms);
    else
        fprintf('Files exist %s, %s, %s, %s, and %s...\n',fname_cluster2_sliding_ips,fname_cluster2_times,fname_cluster2_labels,fname_cluster2_waveforms,fname_cluster2_features);
    end;
end;

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');

end

function clips=extract_clips(Y,times,clip_size)
% Extract clips centered at times

[M,N]=size(Y);
T=clip_size;
C=length(times);

clips=zeros(M,T,C);
tt1=-floor((clip_size+1)/2);
tt2=tt1+clip_size-1;
if (min(times+tt1)<1) error('Invalid time in extract_clips'); end;
if (max(times+tt2)>N) error('Invalid time in extract_clips'); end;
for j=1:C
	clips(:,:,j)=Y(:,times(j)+tt1:times(j)+tt2);
end;

end
