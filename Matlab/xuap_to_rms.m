% xuap_to_rms

res = dir('xuap.*');
u_rms = 0;
v_rms = 0;
w_rms = 0;
u_scale = 0;
totalN = 0;

for i  = 1:length(res)
% f = load('F:\Two-phase\Set_B\11Sep_run3\res\xuap.102711');
    f = load(res(i).name);
    u = f(:,9);
    v = f(:,10);
    w = f(:,11);
    
    u_rms = u_rms + sum(u.^2);
    v_rms = v_rms + sum(v.^2);
    w_rms = w_rms + sum(w.^2);
    totalN = totalN + size(f,1);
end
u_rms = sqrt(u_rms/totalN);
v_rms = sqrt(v_rms/totalN);
w_rms = sqrt(w_rms/totalN);

u_scale = 1/3*sqrt((u_rms + v_rms + w_rms)/(totalN));

% hold on
% hist(f(:,9))
% hist(f(:,10))
% hist(f(:,11))
