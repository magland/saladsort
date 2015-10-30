function [X t m info] = ss_detectevents(A,opts)
% A snapshot of ahb's detectevents - 5/22/2015 - replacing d by A, removing
% dependence on opts.Twin and d.samplefreq
%
% MAJOR CHANGE:::: the output t corresponds to the centers of the events,
% not the beginnings of the clip windows (sorry ahb that this may cause
% trouble!!), except this does not apply for multiple events (not really
% using multiple events at this stage much anymore, I think)
%
% SS_DETECTEVENTS - detect single and multple spiking events on a small # channels
%
% [X t m] = detectevents(A) where A is the MxN data array (filtered)
%  extracts believed single spiking events in s, with timeshifts t, and
%  multiple-spiking events in a cell array m.
%  Small number of channels only, since detects across all channels.
%
%  jfm notes: If opts.samplefreq is not specified, then we default to 20000
%  and we do not attempt to detect multiple events -- because that really
%  depends on dt.
%
% [X t m info] = detectevents(A,opts) controls certain options and outputs info.
%   opts.verb = 0 (silent), 1 diagnostic output (default)
%   opts.Nt = size of window in timepoints (default 30) (jfm)
%     or if opts.Twin and opts.samplefreq are specified,
%     then opts.Nt=ceil(opts.Twin*opts.samplefreq)
%   opts.quan = detection event threshold quantile (default 0.01, ie 1%)
%               (note 0.005 better for 'b', 0.01 for 'e'). Overriden by thresh.
%   opts.thresh = threshold for detection signal, overrides opts.quan.
%   opts.sig = 'a' detection signal = "energy" u^+(tau.u_t)^2
%              'm' detection signal = minimum over channels
%
% Outputs:
% X - 3d array (M*Nt*Ns) of purported single-spike events (same as Jeremy)
%      Nt is number of window samples, same for all events. Ns = # spike events
% t - vector (1*Ns) of indices of starts of single-spike windows (in 1 to N)
% m - struct for "multiple" events:
%     m.t - (1*Nm) start time indices (in samples, 1 to N)
%     m.Ts - (1*Nm) lengths of each window in samples
%     m.X - concatenated M*m.Ttot signal array (m.Ttot = sum(m.Ts))
%     m.tptr - (1*Nm) column pointers to start of each event in m.X
%     m.Ns - stores Nm

% todo: 1) make T1spike opts.
% 2) Break out a detection signal func for detectevents and alignspikes ?
% No point if align doesn't use it.
%
% Barnett 11/17/14
% 12/4/14: tidied wj0,wj1 detection, fixed spill-over end of 1:N samples
% 2/10/15: opts for plain thresh, channelwise-min. Changed multiple m format
% 2/19/15: fixed offset by 1 bug in m.X
% 3/12/15: o.sig default now 'm'

%addpath([fileparts(mfilename('fullpath')),'/internal']);

if nargin<1, test_ss_detectevents; return; end

if nargin<2, opts = struct; end
if ~isfield(opts,'verb'), opts.verb = 1; end
can_detect_multiple_events=1; %jfm
if (isfield(opts,'samplefreq'))
	%only default the Twin if we know the sample frequency (because
	%otherwise we default the samplefreq, which may not be the true
	if (~isfield(opts,'Twin'))
		opts.Twin = 0.003; % width of a single-event window (s), controls Nt
	end;
	if (~isfield(opts,'Nt'))
		opts.Nt=ceil(opts.Twin*opts.samplefreq); %jfm
	end;
else
	can_detect_multiple_events=0; %we can't do this if we don't know samplefreq
	opts.samplefreq=20000; %I don't think this param is really needed, except in computing the sig below, probably doesn't matter
end;
if ~isfield(opts,'quan'), opts.quan = 0.01; end % a threshold parameter: top 1% of data (~ expected overall firing rate?)
if ~isfield(opts,'sig'), opts.sig = 'm'; end

opts.dt=1/opts.samplefreq;
if (~isfield(opts,'Nt'))
	if (isfield(opts,'Twin'))
		opts.Nt=ceil(opts.Twin*opts.samplefreq); %jfm
	else
		opts.Nt=60;
		opts.Twin=opts.dt*opts.Nt; %jfm
	end;
end;

Nt = opts.Nt;  % # samples in single-event window
[M N] = size(A); 

% one-channel signal for detection...
if strcmp(opts.sig,'a')  % Alex's method: sqrt(energy)
  Ap = [diff(A,[],2) zeros(M,1)]; % upwind stencil derivative
  tau = 0.0002;    % spike rise timescale (sec); hopefully dataset-independent
  sig = sqrt(sum(A.^2 + (tau*opts.samplefreq)^2*Ap.^2,1)); % energy: u^2 + (tau.u')^2
  %figure; hist(sig,50); figure; semilogy((1:N)*d.dt,sig)
elseif strcmp(opts.sig,'m')
 sig = abs(min(A,[],1));     % abs of minimum over channels
end

% Threshold this one-channel signal...
ssig = sort(sig);
if ~isfield(opts,'thresh')
  thresh = ssig(round((1-opts.quan)*N)); % use quantile for threshold for sig
else, thresh = opts.thresh; end
if opts.verb, fprintf('threshold = %.3g\n',thresh); end
% todo: auto thresh, eg so many std devs below mean?

% fatten the detection peaks to include Nt/2 either side, choose centering...
b = [0, conv(double(sig>thresh), ones(1,Nt),'full'), 0]; % 0-padded spread peaks
sh = -floor(Nt/2);       % shift, adjust to center window (to nearest sample)
b = diff(b>0); wj0 = find(b==1)+sh; wj1 = find(b==-1)+sh; % inds starts & stops
t = round((wj0+wj1)/2 - Nt/2); % start indices (in samples 1...N) if single
falloff = (t<1 | t+Nt-1>N);  % true for windows who'll fall off even if short
t = t(~falloff);          % kill them
wj0 = wj0(~falloff); wj1 = wj1(~falloff); % (for wj0 could leave <1 or wj1 >N)
info.b=b; info.sig=sig; info.thresh=thresh; info.wj0=wj0; info.wj1=wj1; % diagn
Ne = numel(wj0); % total # events

% use event durations to split into single- vs multiple-spike events
T1spike = 0.001;      % max time a single spike can be above-threshold (sec)
if (can_detect_multiple_events)
	sing = (wj1-wj0)*opts.dt < (opts.Twin+T1spike);     % duration criterion
else
	sing = ones(size(t)); %jfm -- all single
end;
js = find(sing); t = t(js); % "single" inds in window list, their start times
Ns = numel(js);  % # purported single-spike events
Nm = Ne-Ns;    % # multi-spike events (all are kept)
if opts.verb, fprintf('detectevents: Nt=%d, %d single-events & %d multi-events\n',Nt,Ns,Nm); end
info.Nt = Nt;

%Here's the craziness by jfm
t=round(t+Nt/2);

X = nan(M,Nt,Ns);
for n=1:Ns         % loop over single events
  j = t(n) + (0:Nt-1);       % sample indices in window
  X(:,:,n) = A(:,j);       % copy event into output array
end

js = find(~sing); % probably multiple-spike events...
m.t = max(1,wj0(js));
m.tlast = min(N,wj1(js)-1);      % last index of each window
m.Ts = m.tlast - m.t;            % lengths in signals
m.Ttot = sum(m.Ts);
m.Ns = numel(m.Ts);
m.tptr = [1, 1+cumsum(m.Ts(1:Nm-1))]; % pointers to start col indices in m.X
m.X = nan(M,m.Ttot);
for n=1:Nm, j=0:m.Ts(n)-1; m.X(:,m.tptr(n)+j) = A(:,m.t(n)+j); end
%%%%%%%%%

function local_show_detect(d,t,m,info)
% SHOW_DETECT - overlay detection windows (single, multiple) on signal plot
%
% show_detect(d,t,m,info)
% d - raw EC data object
% t - list of window start times as in detectevents
% m - struct containing multiple-spike events
% info - as output by detectevents
%
% Shows single-tagged events w/ pink rectangle, multiple-ones w/ green.

% Barnett 2/10/15. taken from test_detectevents.

Nt = info.Nt;
Ns = numel(t);  % # single events
local_viewraw(d);    % time axis in physical units
ax = get(gcf,'children'); set(gcf,'currentaxes',ax(2)); % goto other axes
hold on; v = axis; ylo = v(3); yhi = v(4); % get y range for vertical bars
patch([t;t+Nt;t+Nt;t]*d.dt, repmat([yhi;yhi;ylo;ylo],[1 Ns]), ...
      [1 .8 .8], 'linestyle','none');    % single-spike

for i=1:numel(m.t), t = [m.t(i) m.tlast(i)];          % multi-spike [start stop]
  patch([t(1);t(2);t(2);t(1)]*d.dt, [yhi;yhi;ylo;ylo], ...
      [.8 1 .8], 'linestyle','none');
end

sc = 0.1*(yhi-ylo)/info.thresh; % also show threshhold criterion
N = size(d.A,2);
hold on; plot((1:N)/d.samplefreq, info.sig*sc,'k-'); hline(info.thresh*sc);

%% For testing purposes!!
function fig = local_viewraw(d,fig)
% VIEWRAW - plot raw data and electrode picture
%
% viewraw(d) where d is raw data struct opens new figure with signals
% viewraw(d,fig) overplots in figure fig
%
% viewraw with no arguments runs a test.
%
% Part of sspack.
%
% To do: make interactive
%
% Barnett 11/14/14. removed addpath 12/17/14

if nargin<1, test_viewraw; return; end

sep = 1.5 * max(abs(d.A(:)));  % vertical separation in signal units
[M N] = size(d.A);
t = (1:N)/d.samplefreq;
j = 1:N; % time indices to show: all of it
if nargin<2, figure; else, figure(fig); clf reset; end
plot(t(j),bsxfun(@plus,d.A(:,j)',sep*(0:M-1))); % signal plot, displaced
xlabel('t (s)'); axis tight;
set(gca,'ytick',(0:M-1)*sep,'yticklabel',num2cellstr(1:M)); % label elec #s
title([d.name sprintf(': M=%d N=%d T=%g',M,N,d.T)]);

axes('position',[.85 .85 .15 .15]); % inset
co = get(gca,'colororder'); Nco = size(co,1); % standard color ordering
for m=1:M, c = co(mod(m-1,Nco)+1,:); % 1x3 color vector matching signal graphs
  plot(d.electrodelocs(1,m),d.electrodelocs(2,m),'.','color',c); hold on;
  text(d.electrodelocs(1,m),d.electrodelocs(2,m),sprintf('%d',m),'color',c);
end
set(gca,'xtick',[],'ytick',[]);
axis equal tight; v=axis; axis(v + 0.5*[-1 1 -1 1]); % pad view a bit

function test_ss_detectevents % makes a plot showing event windows
d = ss_load_test_data('e');
A = ss_freqfilter(d.A,d.samplefreq,300,[]);
[s t m info] = ss_detectevents(A);
local_show_detect(d,t,m,info);
set(gca,'xlim',[0 1]); % only 1st 1 sec for sanity

% Now let's try it with samplefreq so it can try to detect multiple events
% as well
[s t m info] = ss_detectevents(A,struct('samplefreq',20000));
local_show_detect(d,t,m,info);
set(gca,'xlim',[0 1]); % only 1st 1 sec for sanity

