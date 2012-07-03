function [xi,yi,uf,vf] = plot_quiver3_xuap(data,varargin)
[xi,yi] = meshgrid(-0.03:.005:.06,-0.03:.005:0.06);
[uf,vf] = deal(zeros(size(xi)));
figure, subplot(121), hold on, box on, axis([-0.022 0.056 -0.022 0.056]), axis square
warning off
for i = 1:length(data)
    if any(data(i).uf)
        try
        uf = uf + griddata(data(i).xf,data(i).yf,data(i).uf,xi,yi,'nearest');
        vf = vf + griddata(data(i).xf,data(i).yf,data(i).vf,xi,yi,'nearest');
        catch
            continue
        end
    end
    uf1 = data(i).uf;
    vf1 = data(i).vf;
    ind = find(uf1 < prctile(uf1,95) & vf1 < prctile(vf1,95));
    ind = ind(1:5:end);
    % ind = 1:length(uf);
    if ~isempty(ind)
        quiver(data(i).xf(ind),data(i).yf(ind),data(i).uf(ind),data(i).vf(ind),0.4);
    end
end
warning on
uf = uf/length(data);
vf = vf/length(data);
subplot(122), quiver(xi,yi,uf,vf); axis tight, axis([-0.022 0.056 -0.022 0.056]), axis square, box on
%{
figure, hold on
for i = 1:length(data)
%     uf = data(i).uf;
%     vf = data(i).vf;
%     wf = data(i).wf;
%     ind = find(uf < prctile(uf,99) & vf < prctile(vf,99));
    % ind = 1:length(uf);
    if ~isempty(ind)
        quiver3(data(i).xf,data(i).yf,data(i).zf,data(i).uf,data(i).vf,data(i).wf);
    end
end
%}
