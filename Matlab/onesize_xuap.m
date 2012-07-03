function xuap = onesize_xuap(xuap)
fields = fieldnames(xuap);
numParticles = zeros(length(xuap),1);
for i = 1:length(xuap)
    numParticles(i) = length(xuap(i).(fields{1}));
end
maxNumParticles = max(numParticles);

for i = 1:length(xuap)
    if numParticles(i) < maxNumParticles
        for j = 1:length(fields)
            xuap(i).(fields{j})(maxNumParticles) = 0;
            xuap(i).t = xuap(i).t(:)';
        end
    end
end
