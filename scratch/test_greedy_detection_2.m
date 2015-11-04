function test_greedy_detection_2

close all;

opts=configuration_opts;
opts.num_timepoints=1e7;
opts.Nt=120;
channels=3:7;

data=grab_data(opts);
X=data.X(channels,:); %restrict to subset of channels
N=size(X,2); % number of timepoints
opts.tt=(1:opts.Nt)-ceil((opts.Nt+1)/2); % eg, [-60:59], for convenience
opts.aa=opts.Nt+1:N-opts.Nt-1; % eg, [121:N-121], for convenience

opts.times_to_exclude=[];
while 1
    [times,avg_waveform]=detect_a_spike(X,opts);
    opts.times_to_exclude=[opts.times_to_exclude,times];
    if (~isempty(avg_waveform))
        ss_view_waveforms(avg_waveform);
        drawnow;
    end;
end;

end

function [times,avg_waveform]=detect_a_spike(X,opts)
fprintf('Excluding %d times...\n',length(opts.times_to_exclude));
for j=1:length(opts.times_to_exclude)
    X(:,opts.times_to_exclude(j)+opts.tt)=0;
end;
fprintf('Computing sliding sumsqr...\n');
N=size(X,2);
sliding_sumsqr=compute_sliding_sumsqr(X,opts);
num_tries=0;
fprintf('Sorting...\n');
[~,sort_inds]=sort(sum(X(:,opts.aa).^2,1));
sort_inds=opts.aa(sort_inds);
kernel_time=sort_inds(end);
kernel=X(:,kernel_time+opts.tt);
for pass=1:2
    fprintf('*Pass %d*\n',pass);
    kernel_sumsqr=sum(kernel(:).^2);
    fprintf('Computing sliding inner product...\n');
    sliding_inner_product=compute_sliding_inner_product(X,kernel,opts);
    sliding_resid=sliding_sumsqr+kernel_sumsqr-2*sliding_inner_product;
    reduction_scores=sliding_sumsqr-sliding_resid;
    fprintf('Finding local global minima...\n');
    candidate_times=find_local_global_maxima(reduction_scores,opts.Nt,kernel_sumsqr*0.05);
    if (length(candidate_times)<10)
        times=candidate_times;
        avg_waveform=[];
        fprintf('Not enough candidate times (%d).\n',length(candidate_times));
        return;
    end;
    if (pass==1)
        [~,inds00]=sort(reduction_scores(candidate_times));
        candidate_times=candidate_times(inds00(end-9:end)); %take the top 10 for the first pass
    end;
    if (isempty(find(candidate_times==kernel_time)))
        candidate_times=[candidate_times,kernel_time];
    end;
    
    clips=extract_clips(X,candidate_times,opts.Nt);
    if (pass==2)
        random_times=randsample(opts.aa,length(candidate_times));
        %make sure we don't duplicate a time!
        tmp=zeros(1,N);
        tmp(random_times)=1;
        tmp(candidate_times)=0;
        random_times=find(tmp);
        noise_clips=extract_clips(X,random_times,opts.Nt);
        clips=cat(3,clips,noise_clips);
    end;
    fprintf('Event features...\n');
    FF=ss_eventfeatures(clips);
    FF=FF(1:3,:);
    fprintf('ISO-SPLIT...\n');
    labels=isosplit(FF);
    fprintf('...\n');
    kernel_ind=find(candidate_times==kernel_time);
    kernel_label=labels(kernel_ind);
    kernel_cluster_indices=find(labels==kernel_label);
    fprintf('Found %d events in cluster %d...\n',length(kernel_cluster_indices),kernel_label);
    if (pass==1)
        kernel=median(clips(:,:,kernel_cluster_indices),3);
    elseif (pass==2)
        noise_ind=length(candidate_times)+1;
        noise_label=labels(noise_ind);
        if (noise_label==kernel_label)
            fprintf('Noise label is kernel label.\n');
            times=kernel_time;
            avg_waveform=[];
            return;
        end;
        ss_view_clusters(FF,labels);
        kernel_cluster_indices=kernel_cluster_indices(find(kernel_cluster_indices<=length(candidate_times)));
        times=candidate_times(kernel_cluster_indices);
        avg_waveform=median(clips(:,:,kernel_cluster_indices),3);
    end;
end;

end

function ret=compute_sliding_sumsqr(X,opts)
N=size(X,2);
ret=zeros(1,N);
for dt=opts.tt
    fprintf('.');
    ret(opts.aa)=ret(opts.aa)+sum(X(:,opts.aa+dt).^2,1);
end;
fprintf('\n');
end

function ret=compute_sliding_inner_product(X,kernel,opts)
N=size(X,2);
ret=zeros(1,N);
for dt=opts.tt
    fprintf('.');
    kernel_vals=repmat(kernel(:,dt-opts.tt(1)+1),1,length(opts.aa));
    ret(opts.aa)=ret(opts.aa)+sum(X(:,opts.aa+dt).*kernel_vals,1);
end;
fprintf('\n');
end

function inds=find_local_global_maxima(A,Nt,minval)
N=length(A);
inds=find(A>minval);
inds=inds(find((inds>Nt)&(inds<=N-Nt)));
use_it=ones(1,length(inds));
for j=1:length(inds)
    k=j-1;
    while (k>=1)&&(inds(k)>=inds(j)-Nt)
        if (A(inds(k))>A(inds(j)))
            use_it(j)=0;
        else
            use_it(k)=0;
        end;
        k=k-1;
    end;
end;
inds=inds(find(use_it));
end

function clips=extract_clips(Y,times,clip_size)
% Extract clips centered at times

[M,N]=size(Y);
T=clip_size;
C=length(times);

clips=zeros(M,T,C);
tt1=-floor(clip_size/2);
tt2=tt1+clip_size-1;
if (min(times+tt1)<1) error('Invalid time in extract_clips'); end;
if (max(times+tt2)>N) error('Invalid time in extract_clips'); end;
for j=1:C
	clips(:,:,j)=Y(:,times(j)+tt1:times(j)+tt2);
end;

end

function data=grab_data(opts)

raw_mda=opts.raw_mda;

%if (~isempty(whos('global','global_X')))
%    global global_X;
%    data.X=global_X;
%else
    fprintf('Reading %s... ',raw_mda);
    timerA=tic;
    
    if (isinf(opts.num_timepoints))
        data.X=readmda(raw_mda);
    else
        data.X=readmda_data_beginning(raw_mda,opts.num_timepoints);
    end;
%    global global_X;
%    global_X=data.X;
    fprintf('\nElapsed: %g seconds',toc(timerA));
    fprintf('\n');
%end;

end
