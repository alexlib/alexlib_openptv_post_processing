function clean_traj_old(xuap)
% added time to xuap


% plot vs time shows that there are many zeros and non-single particle
% frames
% figure, hold on, for i = 1:length(xuap); for j = 1:length(xuap(i).xf), plot(xuap(i).t/50,xuap(i).xf(j),'.'); end; end;

% remove zeros

tmp = xuap;
flds = fieldnames(xuap);


id = [];
frameid = [];
for i = 1:length(xuap)
    for j = 1:length(xuap(i).xf)
        if (xuap(i).xf(j) == 0 || xuap(i).yf(j) == 0 || xuap(i).zf(j) == 0)
            id = cat(1,id,j);
        end
    end
    if ~isempty(id)
        for k = 1:length(flds)-1
            xuap(i).(flds{k})(id) = [];
        end
        %         xuap(i).yf(id) = [];
        %         xuap(i).zf(id) = [];
        id = [];
    end
    if isempty(xuap(i).xf)
        frameid = cat(1,frameid,i);
    end
end
xuap(frameid) = [];
