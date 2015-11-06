function do_view(ch)

if nargin<1 ch=0; end;

%close all;

addpath('view');
addpath('view/internal');
addpath('util');
addpath('spikespy/matlab');

opts=configuration_opts;

data=grab_data(opts);

if (ch>0)
    view_channel_results(ch,opts,data);
else
    view_waveforms(opts);
    view_events(5e6,opts,data);
end;

% Just view the raw data
%spikespy(data.X);



end

function view_events(num_timepoints,opts,data)
times=readmda(opts.consolidate_times_path);
labels=readmda(opts.consolidate_labels_path);
inds=find(times<=num_timepoints);
times=times(inds);
labels=labels(inds);
spikespy({data.X(:,1:min(end,num_timepoints)),times,labels});    
end

function view_waveforms(opts)
waveforms=readmda(opts.consolidate_waveforms_path);
%ss_view_waveforms(waveforms);
%set(gcf,'position',[100,100,1000,1000]);
view_waveforms_spikespy(waveforms);
end

function view_channel_results(ch,opts,data)

figure; position_figure(gcf,'center',1500,500);
subplot(1,3,1);
view_detected_waveforms(ch,opts);
title(sprintf('Detected waveforms for channel %d',ch));
subplot(1,3,2);
view_detected_clusters(ch,opts,1);
title(sprintf('Detected clusters for channel %d (pos)',ch));
subplot(1,3,3);
view_detected_clusters(ch,opts,-1);
title(sprintf('Detected clusters for channel %d (neg)',ch));

PARAMS.ch=ch; PARAMS.opts=opts; PARAMS.data=data;
h=uicontrol(gcf,'Style','pushbutton','String','View events in spikespy','Position',[20,20,300,20],'Callback',{@cb_view_events_in_spikespy,PARAMS});

end

function cb_view_events_in_spikespy(source,callbackdata,PARAMS)
view_detected_events_spikespy(PARAMS.ch,PARAMS.opts,PARAMS.data);
end

function position_figure(fig,pos,W,H)

screen_size = get(0,'screensize'); %The screen size
screen_width = screen_size(3);
screen_height = screen_size(4);

if (strcmp(pos,'center'))
    x=floor(screen_width/2-W/2);
    y=floor(screen_height/2-H/2);
end;

set(fig,'position',[x,y,W,H]);

end










function data=grab_data(opts)

raw_mda=opts.raw_mda;

if (~isempty(whos('global','global_X')))
    global global_X;
    data.X=global_X;
else
    fprintf('Reading %s... ',raw_mda);
    timerA=tic;
    if ~isempty(opts.timepoints)
        data.X=readmda_data_beginning(raw_mda,max(opts.timepoints));
        data.X=data.X(:,opts.timepoints);
    else
        data.X=readmda(raw_mda);
    end;
    global global_X;
    global_X=data.X;
    fprintf('\nElapsed: %g seconds',toc(timerA));
    fprintf('\n');
end;

end
