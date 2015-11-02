function greedy_detection

opts=configuration_opts;
opts.num_timepoints=2e6;

data=grab_data(opts);

ch=3:5;
X=data.X(ch,:);
N=size(X,2);
Nt=120;
tt=(1:Nt)-ceil((Nt+1)/2);
aa=Nt+1:N-Nt-1;

num_repeats=10;
for rr=1:num_repeats

    okay_to_proceed=0;
    while ~okay_to_proceed
        sliding_sumsqr=zeros(1,N);
        fprintf('Computing sliding sumsqr\n');
        for dt=tt
            sliding_sumsqr(aa)=sliding_sumsqr(aa)+sum(X(:,aa+dt).^2,1);
        end;
        [~,sort_inds]=sort(sum(abs(X).^2,1));
        maxtime=sort_inds(ceil(length(sort_inds)*1));
        %[~,maxtime]=max(sliding_sumsqr);
        kernel=X(:,maxtime+tt);
        ss_view_waveforms(kernel); title('kernel');
        kernel_sumsqr=sum(sum(kernel.^2));
        sliding_ip=zeros(1,N);
        fprintf('Computing sliding ip\n');
        for dt=tt
            kernel_vals=repmat(kernel(:,dt-tt(1)+1),1,length(aa));
            sliding_ip(aa)=sliding_ip(aa)+sum(X(:,aa+dt).*kernel_vals,1);
        end;
        %sum(a-b)^2 = Saa+Sbb-2Sab
        resid_sumsqr=sliding_sumsqr+kernel_sumsqr-2*sliding_ip;
        reduction_scores=sliding_sumsqr-resid_sumsqr;
        reduction_times=find(reduction_scores>0);

        %Now make sure we only use those that are global maxima over the radius Nt
        use_it=ones(1,length(reduction_times));
        for j=1:length(reduction_times)
            if ((reduction_times(j)-Nt)<1) use_it(j)=0; end;
            if ((reduction_times(j)+Nt)>N) use_it(j)=0; end;
            k=j-1;
            while (k>=1)&&(reduction_times(k)>=reduction_times(j)-Nt)
                if (reduction_scores(reduction_times(k))>reduction_scores(reduction_times(j)))
                    use_it(j)=0;
                else
                    use_it(k)=0;
                end;
                k=k-1;
            end
        end;
        reduction_times=reduction_times(find(use_it));

        if (length(reduction_times)>=10)
            okay_to_proceed=1;
        else
            fprintf('Not enough spikes (%d), trying again.\n',length(reduction_times));
            okay_to_proceed=0;
            X(:,maxtime+tt)=0;
        end
    
    end

    reduction_times=[reduction_times,randsample(N,length(reduction_times))'];

    fprintf('Extracting clips...\n');
    clips=extract_clips(X,reduction_times,Nt);

    fprintf('ISO-SPLIT...\n');
    FF=ss_eventfeatures(clips);
    FF=FF(1:3,:);
    labels=isosplit(FF);

    spikespy({X,reduction_times,labels});

    fprintf('Computing average waveforms...\n');
    K=max(labels);
    WF=zeros(size(X,1),Nt,K);
    for k=1:K
        inds=find(labels==k);
        %WF(:,:,k)=mean(clips(:,:,inds),3);
        WF(:,:,k)=median(clips(:,:,inds),3);
    end;

    ind0=find(reduction_times==maxtime);
    label0=labels(ind0);
    WF0=WF(:,:,label0);
    times0=reduction_times(find(labels==label0));

    ss_view_waveforms(WF);
    disp(label0)
    
    for t0=times0
        X(:,t0+tt)=0;
    end;

    ss_view_clusters(FF,labels);
    drawnow;
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
