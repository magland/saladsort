function step2_detect(opts,data)

timerA=tic;

fprintf('Step 2: Detect... ');

detect_times_prefix=opts.detect_times_prefix;
detect_clips_prefix=opts.detect_clips_prefix;

X=data.X;

fprintf('Detecting...\n');
oo.thresh=5;
oo.Nt=opts.clip_size;
for j=1:size(X,1)
    fname=[detect_times_prefix,sprintf('%d.mda',j)];
    fname2=[detect_clips_prefix,sprintf('%d.mda',j)];
    if ((~exist(fname,'file'))||(~exist(fname2,'file')))
        T=detect_events(X(j,:),oo);
        fprintf('Writing %s (%d events detected)...\n',fname,length(T));
        writemda(T,fname);
        
        clips=extract_clips(X,T,oo.Nt);
        fprintf('Writing %s...\n',fname2);
        writemda(clips,fname2);
    else
        fprintf('Files exists %s and %s...\n',fname,fname2);
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
tt1=-ceil((clip_size)/2);
tt2=tt1+clip_size-1;
if (min(times+tt1)<1) error('Invalid time in extract_clips'); end;
if (max(times+tt2)>N) error('Invalid time in extract_clips'); end;
for j=1:C
	clips(:,:,j)=Y(:,times(j)+tt1:times(j)+tt2);
end;

end

function T=detect_events(X,opts)

Nt=opts.Nt;
thresh=opts.thresh;

absX=abs(X);
absX=max(absX,1); %max over channels (but best to put in one channel at a time)
times=find(absX>thresh);

%Now make sure we only use those that are global maxima over the radius
%Nt/2
Nt_2=ceil(Nt/2);
use_it=ones(1,length(times));
for j=1:length(times)
    if ((times(j)-Nt_2)<1) use_it(j)=0; end;
    if ((times(j)+Nt_2)>size(X,2)) use_it(j)=0; end;
    k=j-1;
    while (k>=1)&&(times(k)>=times(j)-Nt_2)
        if (absX(times(k))>absX(times(j)))
            use_it(j)=0;
        else
            use_it(k)=0;
        end;
        k=k-1;
    end
end;

T=times(find(use_it));

end

