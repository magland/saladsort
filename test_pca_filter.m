npca=10;
M=10;
T=80;
NC=100;
clips=rand(M,T,NC);

[FF,info]=ss_eventfeatures(clips);
FF=FF(1:npca,:);
SS=info.subspace;
SS=SS(:,:,1:npca);

SS0=reshape(SS,[M*T,size(SS,3)]);
clips0=reshape(clips,[M*T,NC]);
clips_filt=reshape(SS0*FF,[M,T,NC]);
disp(info);

resid=clips-clips_filt;
sqrt(sum(resid(:).^2))