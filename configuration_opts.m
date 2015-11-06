function opts=configuration_opts

opts.raw_dat='raw/ms11d45.dat';
opts.channels=[37:52,68,69];
opts.raw_mda='raw/ms11d45_groupA_pre.mda';
opts.detect_times_prefix='output/detect/times_';
opts.detect_clips_prefix='output/detect/clips_';

opts.cluster_times_prefix='output/cluster/times_';
opts.cluster_labels_prefix='output/cluster/labels_';
opts.cluster_waveforms_prefix='output/cluster/waveforms_';
opts.cluster_features_prefix='output/cluster/features_';

opts.cluster2_sliding_ips_prefix='output/cluster2/sliding_ips_';
opts.cluster2_times_prefix='output/cluster2/times_';
opts.cluster2_labels_prefix='output/cluster2/labels_';
opts.cluster2_waveforms_prefix='output/cluster2/waveforms_';
opts.cluster2_features_prefix='output/cluster2/features_';

opts.consolidate_times_path='output/consolidate/times.mda';
opts.consolidate_labels_path='output/consolidate/labels.mda';
opts.consolidate_waveforms_path='output/consolidate/waveforms.mda';
opts.consolidate_load_channels_path='output/consolidate/load_channels.mda';

opts.fit_times_prefix='output/fit/times_';
opts.fit_labels_prefix='output/fit/labels_';
opts.locations='raw/locations.mda';
opts.adjacency='raw/adjacency.mda';

%opts.timepoints=1:19e6; %Something seems to change around timepoint 19-20 million
opts.timepoints=[1:2.9e6,3.0e6:19e6]; %Something seems to change around timepoint 19-20 million
%opts.timepoints=1:5e6;
opts.num_cluster_features=6;
opts.clip_size=80;
opts.detection_threshold=3;
opts.detection_interval=80;

end
