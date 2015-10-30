function data=step1_preprocess(opts)

fprintf('Step 1: Preprocess... ');
channels=opts.channels;

timerA=tic;

fname_in=opts.raw_dat;
fname_out=opts.raw_mda;
channels=opts.channels;

if (exist(fname_out,'file'))
    fprintf('File %s already exists.\n',fname_out);
    if (~isempty(whos('global','global_X')))
        global global_X;
        data.X=global_X;
    else
        global global_X;
        fprintf('Reading %s... ',fname_out);
        data.X=readmda(fname_out);
        global_X=data.X;
        fprintf('\nElapsed: %g seconds',toc(timerA));
        fprintf('\n');
    end;
    fprintf('%d chan, %g timepoints\n',size(data.X,1),size(data.X,2));
    return;
end;

fprintf('Creating %s...',fname_out');

N=inf;

fprintf('\nReading %s... ',fname_in);
F=fopen(fname_in,'rb');
X=fread(F,[72,N],'int16');
fclose(F);

%Extract group 1
fprintf('Extracting group 1... ');
X=X(channels,:);

fprintf('%d chan, %g timepoints\n',size(X,1),size(X,2));

fprintf('Filtering... ');
X=ss_freqfilter(X,30000,300,10000);
d.A=X; d.samplefreq=30000; d.dt=1/d.samplefreq;
fprintf('Prewhitening... ');
d=channelprewhiten(d,[]);
X=d.A;

fprintf('Normalizing... ');
for j=1:size(X,1)
    stdev=sqrt(var(X(j,:)));
    X(j,:)=X(j,:)/stdev;
    X(j,:)=min(20,max(-20,X(j,:)));
end;
for j=1:size(X,1)
    stdev=sqrt(var(X(j,:)));
    X(j,:)=X(j,:)/stdev;
end;

fprintf('Writing %s... ',fname_out);
writemda(X,fname_out);

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');

data.X=X;

end
