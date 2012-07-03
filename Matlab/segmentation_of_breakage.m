%% Image segmentation
close all
clear all

first=19726;
last=19731;
name_root=['.\WF\'];
% crop sizes
rect = [505.5100  300.5100  143.9800  192.9800;...
    539.5100  312.5100  208.9800  355.9800;...
    571.5100  314.5100  259.9800  396.9800;...
    421.5100  266.5100  230.9800  409.9800];

stats1 = cell(4,last-first+1);
stats2 = stats1;
verbose = 0; % change to True if you want to see a moving particle
for img = 1:4
    % img=1;
    j = 0;
    for i=first:last
        f1=imread([name_root,'img\cam',int2str(img),'.',int2str(i)]);
        j = j + 1;
        [stats1{img,j},stats2{img,j}] =  colloid_image_segmentation_inpoly(f1,rect(img,:),verbose);
        if verbose
            title(sprintf('image %d, frame %d',img,i))
        end
        


    end
end

%% using InPolygon idea:
cand = {};
for img = 1:4
    breakage_found = 0;
    i = first;
    while ~breakage_found
        i = i + 1;
        % for i = first:last
        old = stats1{img,i-first+1};
        %     now1 = stats1{img,i-first};
        new = stats2{img,i-first+1};
        for k = 1:length(old)
            tmp = old(k).ConvexHull; % polygon
            x1 = old(k).Centroid(1);
            y1 = old(k).Centroid(2);
            for m = 1:length(new)
                x21 = new(m).Centroid(1);
                y21 = new(m).Centroid(2);
                for n = m+1:length(new) % all pairs
                    x22 = new(n).Centroid(1);
                    y22 = new(n).Centroid(2);
                    if inpolygon(x21,y21,tmp(:,1),tmp(:,2)) && inpolygon(x22,y22,tmp(:,1),tmp(:,2))
                        j = j + 1;
                        cand{j} = [x1,y1,x21,y21,x22,y22,i];
                        name=[name_root,'img\cam',num2Str(img),'.',int2str(i)];
                        f1=imread(name);
                        colloid_image_segmentation_inpoly(f1,rect(img,:),1);
                        title(sprintf('image %d, frame %d',img,i))
                        hold on
                        plot(x1,y1,'ro',x21,y21,'yx',x22,y22,'yx');
                        hold off
                        breakage_found = 1;
                    end
                end
            end
        end
    end
end

% figure, hold on
% for i = 1:length(cand), plot(cand{i}(1),cand{i}(2),'r.',cand{i}(3),cand{i}(4),'b.'); end
%%
%{
j = 0;
cand = {};
img = 1;
for i = first+1:last
    old = stats2{img,i-first};
    now1 = stats1{img,i-first};
    new = stats2{img,i-first+1};
    for p = 1:length(now1)
        xp = now1(p).Centroid(1);
        yp = now1(p).Centroid(2);
        for k = 1:length(old)
            x1 = old(k).Centroid(1);
            y1 = old(k).Centroid(2);
            if  pdist([xp,yp;x1,y1]) < 1
                for m = 1:length(new)-1
                    x21 = new(m).Centroid(1);
                    y21 = new(m).Centroid(2);
                    for n = m+1:length(new) % all pairs
                        x22 = new(n).Centroid(1);
                        y22 = new(n).Centroid(2);
                        %                 pdist([x21,y21;x22,y22])

                        if pdist([x21,y21;x22,y22]) < 40 && pdist([x1,y1;mean([x21,x22]),mean([y21,y22])]) < 3
                            j = j + 1;
                            cand{j} = [x1,y1,x21,y21,x22,y22,i];
                            name=[name_root,'img\cam',num2Str(img),'.',int2str(i-1)];
                            f1=imread(name);
                            colloid_image_segmentation(f1,rect(img,:),1);
                            title(sprintf('frame %d',i-1))
                            hold on
                            plot(x1,y1,'ro');
                            hold off

                            name=[name_root,'img\cam',num2Str(img),'.',int2str(i)];
                            f1=imread(name);
                            colloid_image_segmentation(f1,rect(img,:),1);
                            title(sprintf('frame %d',i))
                            hold on
                            plot(x1,y1,'ro');
                            plot([x21,x22],[y21,y22],'b-s');
                            hold off
                        end
                    end
                end
            end
        end
    end
end
%}
%%
%{
%%j = 0;
cand = {};
img = 1;
for i = first+1:last
    old = stats2{img,i-first};
    now1 = stats1{img,i-first};
    new = stats2{img,i-first+1};
    for p = 1:length(now1)
        xp = now1(p).Centroid(1);
        yp = now1(p).Centroid(2);
        for k = 1:length(old)
            x1 = old(k).Centroid(1);
            y1 = old(k).Centroid(2);
            if  pdist([xp,yp;x1,y1]) < 2
                %                 for m = 1:length(new)-1
                %                     x21 = new(m).Centroid(1);
                %                     y21 = new(m).Centroid(2);
                %                     for n = m+1:length(new) % all pairs
                %                         x22 = new(n).Centroid(1);
                %                         y22 = new(n).Centroid(2);
                %                         %                 pdist([x21,y21;x22,y22])
                %
                %                         if pdist([x21,y21;x22,y22]) < 40 && pdist([x1,y1;mean([x21,x22]),mean([y21,y22])]) < 3
                %                             j = j + 1;
                %                             cand{j} = [x1,y1,x21,y21,x22,y22,i];
                name=[name_root,'img\cam',num2Str(img),'.',int2str(i-1)];
                f1=imread(name);
                colloid_image_segmentation(f1,rect(img,:),1);
                title(sprintf('frame %d',i-1))
                hold on
                plot(x1,y1,'ro',xp,yp,'yx');
                hold off

                %                             name=[name_root,'img\cam',num2Str(img),'.',int2str(i)];
                %                             f1=imread(name);
                %                             colloid_image_segmentation(f1,rect(img,:),1);
                %                             title(sprintf('frame %d',i))
                %                             hold on
                %                             plot(x1,y1,'ro');
                %                             plot([x21,x22],[y21,y22],'b-s');
                %                             hold off
                %                         end
            end
        end
    end
end
% end
% end
%}

