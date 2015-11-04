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
opts.cluster_times_path='output/cluster/times.mda';
opts.cluster_labels_path='output/cluster/labels.mda';
opts.cluster_waveforms_path='output/cluster/waveforms.mda';
opts.cluster_load_channels_path='output/cluster/load_channels.mda';
opts.fit_times_prefix='output/fit/times_';
opts.fit_labels_prefix='output/fit/labels_';
opts.locations='raw/locations.mda';
opts.adjacency='raw/adjacency.mda';

opts.timepoints=1:19e6; %Something seems to change around timepoint 19-20 million
opts.num_cluster_features=3;
opts.clip_size=80;

end
