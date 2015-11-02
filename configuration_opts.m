function opts=configuration_opts

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

end
