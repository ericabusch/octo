%% Step 2 in the arm tracking process. 
% Bring in first skeletonization of the arm and remove branches to yield single arm.
% eb version
%% Import first skeletonization video

%this code will improve the results of the skeletonization.

imp = 1; %import new file?

% open videos
if imp == 1
    clear all
    % determine what iteration you are running of this script
    % automatically imports the most recent skeletonization avi and
    % creates next 
    iter_num = input('what iteration of skeleton video are you importing? ');
    disp('select the video you want to clean up.');
    [fn, dn] = uigetfile('*.avi');
    destdir = dn;
    iter_num = iter_num + 1;
    fprintf('running skeletonization iteration number %d.\n', iter_num);
    %132-
    root_dir = cd;
    v = VideoReader(sprintf('%s/%s',destdir, fn));
    p = 0;
    while hasFrame(v)
        p = p + 1;
        video1 = readFrame(v);
        video(:,:,p) = double(video1(:,:,1));   
    end

    video(video < 200) = 0;
    video(video > 0) = 1;
end



%% First, identify base and tip of arm 

%IMPORTANT: The first thing you click on should be a spot near the base of the
%arm. The tip location is not used later in the code as of yet.

clear x y
disp('identify base and tip of arm');
str = 'IMPORTANT: The first thing you click on should ';
str = [str newline 'be a spot near the base of the arm. The tip location'];
str = [str newline 'is not used later in the code as of yet.'];
disp(str);
disp('then press return key.');
    
figure, imshow(video(:,:,1));

[x, y] = getpts;

hold on

plot(x(1),y(1),'b*')
disp('the blue * is the base and the red * is the tips');
plot(x(2),y(2),'r*')

%% run on video

for yb = 1:size(video,3)
    
    if yb == 1
        b1b = video(:,:,yb);
        b1b = logical(b1b);
        
        %fill holes to eliminate loops
        filled = imfill(b1b, 'holes');
        
        % image matrix
        b1b = bwmorph(filled, 'skel', Inf);

        %find branch and end points
        B = bwmorph(b1b, 'branchpoints');
        E = bwmorph(b1b, 'endpoints');
        E1 = E;
        [y1, x1] = find(E1);
        dsb = sqrt((y1-y(1)).^2 + (x1-x(1)).^2);
        y1b = y1(dsb == min(dsb));
        x1b = x1(dsb == min(dsb));
        Db = bwdistgeodesic(b1b, x1b, y1b);
        B_loc = find(B);
        dsb1 = Db(E);
        y2b = y1(dsb1 == max(dsb1));
        x2b = x1(dsb1 == max(dsb1));
        
        figure, imagesc(Db)
        hold on
        plot(x1b, y1b,'c*','linewidth',3)
        plot(x2b, y2b,'r*','linewidth',3) 
        
        if length(x1) > 2
            for c = 1:5
                if c==1
                    Dmask = b1b;
                    b1d = b1b;
                else
                    clear x1 y1
                    Dmask = b1d;
                    E = bwmorph(Dmask, 'endpoints');
                    B = bwmorph(Dmask, 'branchpoints');
                    Db = bwdistgeodesic(Dmask,x1b,y1b);
                    Db(isinf(Db)) = 0;
                    B_loc = find(B);
                    [y1,x1] = find(E);
                end

                if length(x1) > 2
                    for k = 1:numel(x1)
                        D = bwdistgeodesic(logical(b1d),x1(k),y1(k));
                        distanceToBranchPt = min(D(B_loc));
                        dst_to_base = Db(y1(k),x1(k));

                        if length(distanceToBranchPt)<1
                            distanceToBranchPt = 1000;
                        end
                        
                        Db(isinf(Db)) = 0;
                        D(isinf(D)) = 0;

                        if dst_to_base<max(max(Db)) && dst_to_base>distanceToBranchPt
                        Dmask(D < distanceToBranchPt) = false;
                        Dmask = bwmorph(Dmask, 'bridge');
                        end
                    end
                    
                    b1d = Dmask;
                    b1d = bwmorph(Dmask, 'bridge');
                    % filled = imfill(b1d, 'holes');
                    % b1d = bwmorph(filled,'skel',Inf);
                    imagesc(b1d,[0 1])
                    drawnow
                    pause(0.03)
                    %Et = bwmorph(b1d, 'endpoints');
                    %[y1t,x1t] = find(Et);  
                end
             %b1d = bwmorph(b1d,'skel',Inf);
            end
        else
            b1d = b1b;
        end  
        % B1D(:,:,yb) = b1d;
    else        
        b1b = video(:,:,yb);
        b1b = logical(b1b);
        %fill holes to eliminate loops
        filled = imfill(b1b, 'holes');
        b1b = bwmorph(filled,'skel',Inf);
        %find branch and end points
        B = bwmorph(b1b, 'branchpoints');
        E = bwmorph(b1b, 'endpoints');
        [y1,x1] = find(E);
        dsb = sqrt((y1-y(1)).^2+(x1-x(1)).^2);
        y1b = y1(dsb==min(dsb));
        x1b = x1(dsb==min(dsb));
        Db = bwdistgeodesic(b1b,x1b,y1b);   
        B_loc = find(B);
        dsb1 = Db(E);         
        y2b = y1(dsb1==max(dsb1));
        x2b = x1(dsb1==max(dsb1));
        
        if length(x1) > 2
            for c = 1:5
                if c==1
                    Dmask = b1b;
                    b1d = b1b;
                else
                    clear x1 y1
                    Dmask = b1d;
                    E = bwmorph(Dmask, 'endpoints');
                    B = bwmorph(Dmask, 'branchpoints');
                    B_loc = find(B);
                    [y1,x1] = find(E);
                    dsb = sqrt((y1-y(1)).^2+(x1-x(1)).^2);
                    y1b = y1(dsb==min(dsb));
                    x1b = x1(dsb==min(dsb));
                    Db = bwdistgeodesic(Dmask,x1b,y1b);
                end
                if length(x1) > 2
                    for k = 1:numel(x1)
                        D = bwdistgeodesic(logical(b1d),x1(k),y1(k));
                        distanceToBranchPt = min(D(B_loc));
                        dst_to_base = Db(y1(k),x1(k));
                        Db(isinf(Db)) = 0;
                        D(isinf(D)) = 0;
                        if dst_to_base<max(max(Db)) && ... 
                                dst_to_base > distanceToBranchPt
                            Dmask(D < distanceToBranchPt) = false;
                            Dmask = bwmorph(Dmask, 'bridge');
                        end
                    end
                    b1d = bwmorph(Dmask, 'bridge');
                end
            end
        else
            b1d = b1b;
        end  
    end
    imagesc(b1d+b1b,[0 2])
    drawnow
    pause(0.03) 
    b1d = bwmorph(b1d, 'bridge');
    B1D(:,:,yb) = b1d;    
end

%% save the output of the reduced skeleton procedure
%NB: recommend to repeat this m-file at least once to generate a
%skel2, then run code on skel2 to get skel3 AVI output. Skel3 is generally close to perfect.

clear m
vid_name = sprintf('%s/skeleton%d.avi', destdir, iter_num);
VidObj = VideoWriter(vid_name);
fprintf('making video %s\n', vid_name);
VidObj.FrameRate = 30;
 
 figure,
 set(gca,'Color','k')
 set(gcf,'Color','k')
open(VidObj);

for b = 1:size(B1D, 3)
    imshow(B1D(:,:,b));
    writeVideo(VidObj, double(B1D(:,:,b)));
end
close(VidObj);
clear VidObj

mat_name = sprintf('%s/skeleton%d.mat', destdir, iter_num);
fprintf('making matfile %s\n', mat_name);
save(mat_name, 'B1D')





