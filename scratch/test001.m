function test001

addpath('view');
addpath('view/internal');
addpath('util');
addpath('spikespy/matlab');

opts=configuration_opts;

data=grab_data(opts);

WW=readmda(opts.cluster_waveforms_path);
[M,T,K]=size(WW);

WW_sparse=zeros(size(WW));
maxval=max(WW(:));
for k=1:K
    W0=WW(:,:,k);
    W0=W0.*(abs(W0)>=maxval*0.05);
    WW_sparse(:,:,k)=W0;
end;

%view_waveforms_spikespy(WW);

[FF,finfo]=ss_eventfeatures(WW);
FF=FF(1:3,:);
subspace=finfo.subspace;
ss_view_clusters(FF,ones(1,size(FF,2)));

return;
%view_waveforms_spikespy(WW_sparse);

WW0=WW_sparse;

AA=zeros(K,K);
for k1=1:K
    W1=WW0(:,:,k1);
    W1=W1/sqrt(W1(:)'*W1(:));
    for k2=1:K
        W2=WW0(:,:,k2);
        W2=W2/sqrt(W2(:)'*W2(:));
        ip=W1(:)'*W2(:);
        AA(k1,k2)=ip;
    end;
end;
%figure; imagesc(AA);

X=data.X;
N=size(X,2);
aa=T:N-T;
t0=ceil((T+1)/2);
Y=zeros(10,N);
parfor k=1:10
    fprintf('.');
    W0=WW_sparse(:,:,k);
    tmp=zeros(1,N);
    for m=1:M
        fprintf('m');
        inds=find(W0(m,:)~=0);
        for ii=1:length(inds)
            tmp(aa)=tmp(aa)+W0(m,inds(ii))*X(m,aa-t0+inds(ii));
        end;
    end;
    tmp=tmp/sqrt(var(tmp));
    Y(k,:)=tmp;
end;
fprintf('\n');

times=readmda(opts.cluster_times_path);
labels=readmda(opts.cluster_labels_path);
spikespy({Y,times,labels});

end

function data=grab_data(opts)

raw_mda=opts.raw_mda;

if (~isempty(whos('global','global_X')))
    global global_X;
    data.X=global_X;
else
    fprintf('Reading %s... ',raw_mda);
    timerA=tic;
    if ~isempty(opts.timepoints)
        data.X=readmda_data_beginning(raw_mda,max(opts.timepoints));
        data.X=data.X(:,opts.timepoints);
    else
        data.X=readmda(raw_mda);
    end;
    global global_X;
    global_X=data.X;
    fprintf('\nElapsed: %g seconds',toc(timerA));
    fprintf('\n');
end;

end
