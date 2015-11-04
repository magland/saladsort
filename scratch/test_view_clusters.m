function test_view_clusters

addpath('view');
addpath('view/internal');
addpath('util');
addpath('spikespy/matlab');

opts=configuration_opts;

data=grab_data(opts);

X=data.X;
M=size(X,1);
N=size(X,2);

times=readmda(opts.consolidate_times_path);
labels=readmda(opts.consolidate_labels_path);
LC=readmda(opts.consolidate_load_channels_path);

k=5; k2=18;
ch=LC(k);
inds=find(labels==k);
times0=times(inds);
clips=extract_clips(X,times0,opts.clip_size);
W=median(clips,3);
NC=size(clips,3);

%inds_noise=randsample(N-opts.clip_size*2,min(N-opts.clip_size*2,NC*3))+opts.clip_size;
inds_noise=find_local_global_maxima(abs(X(ch,:)),80,0);
%inds_noise=find(labels==k2); inds_noise=times(inds_noise);
clips_noise=extract_clips(X,inds_noise,opts.clip_size);
W_noise=median(clips_noise,3);
NC_noise=size(clips_noise,3);

%spikespy(W);
%spikespy(clips);

weight=zeros(size(W));
weight(:,2:end-1,:)=abs(W(:,1:end-2,:))+abs(W(:,2:end-1,:))+abs(W(:,3:end,:));
Wrep=repmat(weight,1,1,NC_noise);

clipsW=clips_noise.*abs(Wrep);
%clipsW=clips;
FFF=ss_eventfeatures(clipsW);
FFF=FFF(1:3,:);
ss_view_clusters(FFF,ones(1,NC_noise));
return;


tmp=(clips-Wrep).*abs(Wrep);
%tmp=(clips-Wrep);
%tmp=squeeze(sqrt(sum(sum(tmp.^2,1),2)));
tmp=squeeze(sum(sum(abs(tmp),1),2));
tmp2=squeeze(sum(sum(clips.*Wrep,1),2))./squeeze(sum(sum(Wrep.*Wrep,1),2));

Wrep=repmat(W,1,1,NC_noise);
tmp_control=(clips_noise-Wrep).*abs(Wrep);
%tmp_noise=(clips_noise-Wrep);
tmp_control=squeeze(sum(sum(abs(tmp_control),1),2));
tmp2_control=squeeze(sum(sum(clips_noise.*Wrep,1),2))./squeeze(sum(sum(Wrep.*Wrep,1),2));

tmp_median=median(tmp);
tmp_rstdev=sqrt(median((tmp-tmp_median).^2));

tmp_control_median=median(tmp_control);
tmp_control_rstdev=sqrt(median((tmp_control-tmp_control_median).^2));

%cutpoint=(tmp_median*tmp_control_rstdev+tmp_control_median*tmp_rstdev)/(tmp_control_rstdev+tmp_rstdev);

alpha=1e10;
mu0=tmp_median; sigma0=tmp_rstdev;
m0=tmp_control_median; s0=tmp_control_rstdev;
H0=log(alpha/sigma0)-log(1/s0); %the a is fixed at 1
a0=1/(2*sigma0^2) - 1/(2*s0^2);
b0=-mu0/(sigma0^2)+m0/(s0^2);
c0=mu0^2/(2*sigma0^2)-m0^2/(2*s0^2)-H0;
discr=sqrt(b0^2-4*a0*c0);
solution1=(-b0-discr)/(2*a0);
solution2=(-b0+discr)/(2*a0);
if ((mu0<=solution1)&&(solution1<=m0))
    cutpoint=solution1;
else
    cutpoint=solution2;
end;
%disp([sigma0,s0,cutpoint]);

%inds_ok=find((tmp-tmp_mean)/tmp_stdev<1.5);
%inds_ok=find_non_outliers(tmp);

figure;
plot(tmp_control,tmp2_control,'k.',tmp,tmp2,'r.');

[~,sort_inds]=sort(tmp);
%spikespy(cat(3,W,clips(:,:,sort_inds(find(tmp(sort_inds)>60)))));
spikespy(cat(3,W,clips(:,:,sort_inds)));

% inds0=find((tmp_control>650)&(tmp2_control>0.7));
% inds0
% spikespy(cat(3,W,clips_noise(:,:,inds0)));

return;

figure;
[nb,xb]=hist(tmp_control,1000);
bh=bar(xb,nb); set(bh,'facecolor','k','edgecolor','k'); hold on;
[nb,xb]=hist(tmp,1000);
bh=bar(xb,nb); set(bh,'facecolor','b','edgecolor','b'); hold on;
%[nb,xb]=hist(tmp(inds_ok),xb);
%bh=bar(xb,nb); set(bh,'facecolor','b','edgecolor','b'); hold off;
h=vline(tmp_median-tmp_rstdev); set(h,'color','b','linewidth',3);
h=vline(tmp_median+tmp_rstdev); set(h,'color','b','linewidth',3);
%h=vline(tmp_control_median-tmp_control_rstdev); set(h,'color','k','linewidth',3);
%h=vline(tmp_control_median+tmp_control_rstdev); set(h,'color','k','linewidth',3);
%h=vline(cutpoint); set(h,'color','r','linewidth',3);
title(sprintf('k=%d, ch=%d',k,ch));

figure;
[nb,xb]=hist(tmp2_control,1000);
bh=bar(xb,nb); set(bh,'facecolor','k','edgecolor','k'); hold on;
[nb,xb]=hist(tmp2,1000);
bh=bar(xb,nb); set(bh,'facecolor','b','edgecolor','b'); hold on;
title(sprintf('k=%d, ch=%d (matched filter)',k,ch));

[~,sort_inds]=sort(tmp2);
%spikespy(cat(3,W,clips(:,:,sort_inds(find(tmp(sort_inds)>60)))));
spikespy(cat(3,W,clips(:,:,sort_inds)));

%[~,sort_inds]=sort(tmp2_control);
%spikespy(cat(3,W,clips_noise(:,:,sort_inds)));

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




% function ret=find_non_outliers(X)
% Xmin=min(X);
% Xmean=mean(X);
% cutoff1=Xmean+(Xmean-Xmin);
% end

function clips=extract_clips(Y,times,clip_size)
% Extract clips centered at times

[M,N]=size(Y);
T=clip_size;
C=length(times);

clips=zeros(M,T,C);
tt1=-ceil((clip_size)/2);
tt2=tt1+clip_size-1;
if (min(times+tt1)<1) error('Invalid time in extract_clips'); end;
if (max(times+tt2)>N) error('Invalid time in extract_clips'); end;
for j=1:C
	clips(:,:,j)=Y(:,times(j)+tt1:times(j)+tt2);
end;

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
