function do_processing

%close all;

%Add some necessary paths
addpath('processing');
addpath('util');
addpath('isosplit');

%Create some directories
if (~exist('output','dir')) mkdir('output'); end;
if (~exist('output/detect','dir')) mkdir('output/detect'); end;

%Set some options and specify the input/output file names
opts=configuration_opts;
opts.verbose=1;

%Preprocess
data=step1_preprocess(opts);
step1a_prepare(opts);

%Detect
step2_detect(opts,data);

%Cluster
step3_cluster(opts,data);

%Final fitting stage needed

end

