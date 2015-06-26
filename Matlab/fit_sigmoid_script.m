% fit sigmoid to buoyancy profiles

% test case
x = z;
figure
dlta = zeros(6,25);
for n = 1:6
    for i = 1:25
        y = B_xmean{n,i};
        
        % figure, plot(x,y);
        
        % scale properly
        mx = mean(x);
        x1 = x - mx;
        
        miny = min(y);
        maxy = max(y);
        y1 = (y - miny)/(maxy-miny);
        
        beta = nlinfit(x1,y1,@sigmoidfit,[1,1,1,1]);
        yhat1 = sigmoidfit(beta,x1);
        [~,k] = min(abs(yhat1 - 0.75));
        
        yhat = yhat1*(maxy-miny) + miny;
        plot(x,y,'b--',x,yhat,'r-',x(k),yhat(k),'ks');
        drawnow
        pause(.5);
        dlta(n,i) = x(k);
    end
end

figure
plot(t,-dlta,'o');


