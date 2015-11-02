function do_view

%close all;

addpath('view');
addpath('view/internal');
addpath('util');
addpath('spikespy/matlab');

opts=configuration_opts;

data=grab_data(opts);

%view_channel_results(3,opts,data);
%view_channel_results(8,opts,data);
view_channel_results(11,opts,data);
%view_channel_results(15,opts,data);

% Just view the raw data
%spikespy(data.X);

end

function view_channel_results(ch,opts,data)

figure; position_figure(gcf,'center',1500,500);
subplot(1,2,1);
view_detected_waveforms(ch,opts);
title(sprintf('Detected waveforms for channel %d',ch));
subplot(1,2,2);
view_detected_clusters(ch,opts);
title(sprintf('Detected clusters for channel %d',ch));
title(sprintf('Channel %d',ch));

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
    data.X=readmda(raw_mda);
    global global_X;
    global_X=data.X;
    fprintf('\nElapsed: %g seconds',toc(timerA));
    fprintf('\n');
end;

end
