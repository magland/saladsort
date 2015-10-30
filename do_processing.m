function do_processing

%close all;

addpath('processing');
addpath('util');
addpath('isosplit');

if (~exist('output','dir')) mkdir('output'); end;
if (~exist('output/detect','dir')) mkdir('output/detect'); end;

opts.verbose=1;
opts.raw_dat='raw/ms11d45.dat';
opts.channels=[37:52,68,69];
opts.raw_mda='raw/ms11d45_groupA_pre.mda';
opts.detect_times_prefix='output/detect/times_';
opts.detect_clips_prefix='output/detect/clips_';
opts.detect_labels_prefix='output/detect/labels_';
opts.detect_waveforms_prefix='output/detect/waveforms_';
opts.detect_features_prefix='output/detect/features_';
opts.locations='raw/locations.mda';
opts.adjacency='raw/adjacency.mda';
opts.num_cluster_features=3;

data=step1_preprocess(opts);
step1a_prepare(opts);
step2_detect(opts,data);
step3_cluster(opts,data);

end

