function view_detected_clusters(ch,opts)

fname=sprintf('%s%d.mda',opts.detect_features_prefix,ch);
try
    FF=readmda(fname);
catch
    error(sprintf('Unable to read file %s. Perhaps you need to run the processing.',fname));
end
fname=sprintf('%s%d.mda',opts.detect_labels_prefix,ch);
labels=readmda(fname);
opts.create_figure=0;
ss_view_clusters(FF,labels,opts);

end