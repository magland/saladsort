function do_processing

%close all;

%Add some necessary paths
addpath('processing');
addpath('util');
addpath('isosplit');

%Create some directories
if (~exist('output','dir')) mkdir('output'); end;
if (~exist('output/detect','dir')) mkdir('output/detect'); end;
if (~exist('output/cluster','dir')) mkdir('output/cluster'); end;
if (~exist('output/cluster2','dir')) mkdir('output/cluster2'); end;
if (~exist('output/consolidate','dir')) mkdir('output/consolidate'); end;
%if (~exist('output/fit','dir')) mkdir('output/fit'); end;

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

%Cluster
%step4_cluster2(opts,data);

step5_consolidate_clusters(opts);

%Consolidate clusters
%step3a_consolidate_clusters(opts,data);

%Final fitting stage needed
%step4_fit(opts,data);

end

