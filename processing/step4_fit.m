function step4_fit(opts,data)

timerA=tic;

fprintf('Step 4: Fit... ');

cluster_waveforms_prefix=opts.cluster_waveforms_prefix;
fit_times_prefix=opts.fit_times_prefix;
fit_labels_prefix=opts.fit_labels_prefix;

X=data.X;

AM=readmda(opts.adjacency);

fprintf('Fitting...\n');
for j=1:size(X,1)
    fname_detect_waveforms=[detect_waveforms_prefix,sprintf('%d.mda',j)];
    fname_fit_times=[fit_times_prefix,sprintf('%d.mda',j)];
    fname_fit_labels=[fit_labels_prefix,sprintf('%d.mda',j)];
    if ((~exist(fname_fit_times,'file'))||(~exist(fname_fit_labels,'file')))
        fprintf('Reading %s... ',fname_detect_waveforms);
        waveforms=readmda(fname_detect_waveforms);
        adjacent_channels=find(AM(:,j));
        XX=X(adjacent_channels,1:1e5);
        WW=waveforms(adjacent_channels,:,:);
        [M,T,K]=size(WW);
        N=size(XX,2);
        fprintf('Computing sliding sumsqr...\n');
        sliding_sumsqr=compute_sliding_sumsqr(XX,opts.clip_size)';
        sliding_inner_products=zeros(N,K);
        sliding_resids=zeros(N,K);
        reduction_scores=zeros(N,K);
        for k=1:K
            fprintf('ip %d/%d...\n',k,K);
            kernel=WW(:,:,k);
            kernel_sumsqr=sum(kernel(:).^2);
            sliding_inner_products(:,k)=compute_sliding_inner_product(XX,kernel);
            sliding_resids(:,k)=sliding_sumsqr+kernel_sumsqr-2*sliding_inner_products(:,k);
            scores=sliding_sumsqr-sliding_resids(:,k);
            scores=scores.*(scores>kernel_sumsqr*0.5);
            reduction_scores(:,k)=scores;
        end;
        [max_reduction_scores,best_labels]=max(reduction_scores,[],2);
        times=find_local_global_maxima(max_reduction_scores,T,0);
        labels=best_labels(times);
        fprintf('Found %d events: ',length(times));
        for k=1:K
            fprintf('%d  ',length(find(labels==k)));
        end;
        fprintf('Writing %s...',fname_fit_times);
        writemda(times,fname_fit_times);
        fprintf('Writing %s...',fname_fit_labels);
        writemda(labels,fname_fit_labels);
    else
        fprintf('Files exist %s and %s...\n',fname_fit_times,fname_fit_labels);
    end;
end;

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');

end

function ret=compute_sliding_sumsqr(X,Nt)
N=size(X,2);
tt=(1:Nt)-ceil((Nt+1)/2); % eg, [-60:59]
aa=Nt+1:N-Nt-1; % eg, [121:N-121]
ret=zeros(1,N);
for dt=tt
    fprintf('.');
    ret(aa)=ret(aa)+sum(X(:,aa+dt).^2,1);
end;
fprintf('\n');
end

function ret=compute_sliding_inner_product(X,kernel)
N=size(X,2);
Nt=size(kernel,2);
tt=(1:Nt)-ceil((Nt+1)/2); % eg, [-60:59]
aa=Nt+1:N-Nt-1; % eg, [121:N-121]
ret=zeros(1,N);
for dt=tt
    fprintf('.');
    kernel_vals=repmat(kernel(:,dt-tt(1)+1),1,length(aa));
    ret(aa)=ret(aa)+sum(X(:,aa+dt).*kernel_vals,1);
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