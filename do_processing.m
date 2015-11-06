function do_processing

%close all;

total_timer=tic;

%Add some necessary paths
addpath('processing');
addpath('util');
addpath('isosplit');

%Create some directories
if (~exist('output','dir')) mkdir('output'); end;
if (~exist('output/detect','dir')) mkdir('output/detect'); end;
if (~exist('output/cluster','dir')) mkdir('output/cluster'); end;
if (~exist('output/consolidate','dir')) mkdir('output/consolidate'); end;
if (~exist('output/residual','dir')) mkdir('output/residual'); end;
%if (~exist('output/fit','dir')) mkdir('output/fit'); end;

%Set some options and specify the input/output file names
opts=configuration_opts;
opts.verbose=1;

%Preprocess
data=step1_preprocess(opts);
fprintf('Bandpass filter...\n');
step1a_prepare(opts);

%Detect
step2_detect(opts,data);

%Cluster
step3_cluster(opts,data);

step4_consolidate_clusters(opts);

step5_residual(opts,data);

fprintf('Total processing time: %.0f sec\n',toc(total_timer));

end

function step5_residual(opts,data)

timerA=tic;

fprintf('Step 5: Residual...\n');

consolidate_times_path=opts.consolidate_times_path;
consolidate_labels_path=opts.consolidate_labels_path;
consolidate_waveforms_path=opts.consolidate_waveforms_path;
consolidate_load_channels_path=opts.consolidate_load_channels_path;
residual_data_path=opts.residual_data_path;
residual_times_path=opts.residual_times_path;
residual_labels_path=opts.residual_labels_path;

AM=readmda(opts.adjacency);
M=size(AM,1);

times=readmda(consolidate_times_path);
labels=readmda(consolidate_labels_path);
waveforms=readmda(consolidate_waveforms_path);
load_channels=readmda(consolidate_load_channels_path);

fprintf('Forming residual...\n');
T=size(waveforms,2);
tt1=-ceil((T)/2);
tt2=tt1+T-1;
tt=tt1:tt2;
X=data.X;
for j=1:length(times)
    X(:,times(j)+tt)=X(:,times(j)+tt)-waveforms(:,:,labels(j));
end;

fprintf('Writing %s...\n',residual_data_path);
writemda(X,residual_data_path);

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');


end


