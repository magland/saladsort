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
        fprintf('Reading %s... ',fname_out);
        if ~isempty(opts.timepoints)
            data.X=readmda_data_beginning(fname_out,max(opts.timepoints));
            data.X=data.X(:,opts.timepoints);
        else
            data.X=readmda(fname_out);
        end;
        global global_X;
        global_X=data.X;
        fprintf('\nElapsed: %g seconds',toc(timerA));
        fprintf('\n');
    end;
    fprintf('%d chan, %g timepoints\n',size(data.X,1),size(data.X,2));
    return;
end;

in_file_ok=1;
if (~exist(fname_in,'file'))
    in_file_ok=0;
else
    dd=dir(fname_in);
    if (dd.bytes~=7800266736)
        fprintf('Wrong number of bytes: %d\n',dd.bytes);
        in_file_ok=0;
    end;
end;

if (~in_file_ok)
    url='http://voms.simonsfoundation.org:50013/rXcGof5b0pDR4zkEDCevBLfo95RB7/frank_ucsf_example1/ms11d45.dat';
    
    fprintf('\n');
    fprintf('It appears that you do not have the raw data file on your computer,\n');
    fprintf('or the file appears to be corrupted or incomplete.\n');
    fprintf('I''ve tried to make it easy for you to obtain it.\n\n');
    
    fprintf('Download the file from here: %s\n',url);
    fprintf('and then place it in the "raw/" directory.\n\n');
    
    fprintf('You could do the following if you had curl:\n');
    fprintf('> curl %s > raw/ms11d45.dat\n\n',url);
    
    fprintf('Once you have obtained the file, re-run the processing script.\n\n');
    
    error('Raw file does not exist. See message above.');
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
