taubins = 0:3:300;
duda = 0*taubins;
count = duda;
theta = duda;


thetabins = -1.1:.1:1.1;
counttheta = 0*thetabins;

for k = 1:length(newtraj)
    traj = newtraj(k);
    for i = 1:length(traj.uf)-1
        for j = i:length(traj.uf)
            u = traj.uf;
            v = traj.vf;
            w = traj.wf;
            ax = traj.axf;
            ay = traj.ayf;
            az = traj.azf;
            tmp = (u(j)-u(i))*(ax(j)-ax(i))+...
                (v(j)-v(i))*(ay(j)-ay(i))+...
                (w(j)-w(i))*(az(j)-az(i));
            tau = j - i;
            ind = bindex(tau,taubins);
            duda(ind) = duda(ind) + tmp;
            count(ind) = count(ind) + 1;
            tmp = cosine([u(j)-u(i),v(j)-v(i),w(j)-w(i)],...
                [ax(j)-ax(i),ay(j)-ay(i),az(j)-az(i)]);
            theta(ind) = theta(ind) + tmp;
            if tau > 5 & tau < 10
                ind = bindex(tmp,thetabins,1);
                counttheta(ind) = counttheta(ind)+1;
            end
        end
    end
end
figure, plot(taubins/50,duda./count);
xlabel('$\tau$, [sec]','Interpreter','Latex');
ylabel('$\langle \Delta u \cdot \Delta a \rangle (\tau) $','Interpreter','Latex')

figure, plot(taubins/50,theta./count);
xlabel('$\tau$, [sec]','Interpreter','Latex');
ylabel('$\langle \theta \rangle (\tau) $','Interpreter','Latex')


figure, plot(thetabins,counttheta./max(counttheta));
