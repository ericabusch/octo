%% Step 3 of the arm morphology tracking process
% erica version
% takes as input the raw video and the final skeletonized arm video

%% Set up: getting inputs, storing raw video in 'video' and skeletonized
% video in 'SKL'
imp = 1; % are you importing a new file?

if imp == 1
    clear all
    sbs1 = 2; %spatial sub-sampling - generally set to 2
    make_changes = 0; % are we making come changes
    % Import raw video file
    disp('Select raw video file.');
    [fn, dn] = uigetfile('*.avi');
    cd(dn)
    
    %132-
    root_dir = cd; %#ok<NASGU>
    v = VideoReader(fn);
    p = 0;
    while hasFrame(v)
        p = p+1;
        video1 = readFrame(v);
        video(:,:,p) = double(video1(1:sbs1:end,1:sbs1:end,1)); %#ok<SAGROW>
    end
    if make_changes == 1
       video = video(:,:,900:1200);
    end
        
    % Import final skeletonized video file
    disp('Select final skeletonized video file.');
    [fn, dn] = uigetfile('*.avi');
    skeldir = dn;
    v = VideoReader(sprintf('%s/%s',skeldir, fn));
    p = 0;
    while hasFrame(v)
        p = p+1;
        video1 = readFrame(v);
        SKL(:,:,p) = double(video1(:,:,1));
    end
    
    SKL(SKL < 200) = 0;
    SKL(SKL > 0) = 1;

end
save_stuff = 1; % do you want to save things as you go?

if save_stuff == 1
    c = date;
    if make_changes == 1
        destdirname = sprintf('NEWocto_arm_tracking_output_%s/',c);
    else
        destdirname = sprintf('octo_arm_tracking_output_%s/',c);
    end
    mkdir(destdirname);
end

%% eliminate spurious waves in tracked data

sr = strel('disk',6);
xb = 0;
for b = 1:size(SKL,3) 
    if b > 1
        skt = sum(sum(video(:,:,b)-video(:,:,b-1)));
    else
        skt = 1;
    end
    
    if skt ~= 0
        xb = xb + 1;
        us1(xb) = b; %#ok<SAGROW>
        for c = 1:3
            if c == 1
                sk = SKL(:,:,b);
            end
            sk = bwmorph(imdilate(sk,sr),'thin','inf');
        end
        imshow(sk);
        drawnow
        SKL(:,:,b) = sk;
        pause(0.02)
    end    
end    
SKL = SKL(:, :, us1);
video = video(:, :, us1);

%% play back raw video with skeletonized arm
clear FRM VidObj

figure,
for b = 1:size(SKL, 3)
    v1 = video(:,:,b)./max(max(video(:,:,b)));
    v2 = SKL(:,:,b);
    frm(:,:,1) = v1 + v2;
    frm(:,:,2) = v1;
    frm(:,:,3) = v1;
    imshow(frm);
    drawnow
    pause(0.02)
    FRM(:,:,:,b) = frm;
    
end

%% select length of arm for analysis
% select lower and higher points that will be displayed layer

%IMPORTANT: You will need to click on three points. The first is the base
%of the arm, the second the farthest segment to include in analysis and the
%third is the midway point, or point that makes contact with the target.

% print to console instructions.
disp('Selecting maximum length for analysis.');
str = 'Please select three points:';
str = [str newline '1. Base of the arm.'];
str = [str newline '2. Farthest segment to include in analysis.'];
str = [str newline '3. Midway point, or point that contacts the target.'];
disp(str);
disp('Then, press the return key.');

SKL = logical(SKL);
% show just the first frame from which to select points
figure, imshow(FRM(:,:,:,1))
[x,y] = getpts;
hold on
plot(x(1), y(1), 'b*')
plot(x(2), y(2), 'r*')
plot(x(3), y(3), 'y*')

% do some math
E = bwmorph(SKL(:, :, 1), 'endpoints');
[y1, x1] = find(E);
[y2, x2] = find(SKL(:, :, 1));
dsb = sqrt((y1 - y(1)).^2 + (x1 - x(1)).^2);
y1b = y1(dsb == min(dsb));
x1b = x1(dsb == min(dsb));
Db = bwdistgeodesic(SKL(:, :, 1), x1b, y1b);
[y2, x2] = find(SKL(:,:,1));
dsb = sqrt((y2 - y(2)).^2+(x2 - x(2)).^2);
y2b = y2(dsb == min(dsb));
x2b = x2(dsb == min(dsb));
y2b = y2b(1);
x2b = x2b(1);
vl = Db(y2b, x2b);
dsb = sqrt((y2 - y(3)).^2 + (x2 - x(3)).^2);
y3b = y2(dsb == min(dsb));
x3b = x2(dsb == min(dsb));
y3b = y3b(1);
x3b = x3b(1);
ml = Db(y3b, x3b);

%% play back the raw video and segmented arm with midpoint
% also save again
clear FRM VidObj

if save_stuff == 1
    
    VidObj = VideoWriter(sprintf('%sraw_with_segmented_arm_30fps.avi', destdirname));
    VidObj = VideoWriter(sprintf('%sraw_segmented_arm_midpoint_30fps.avi', destdirname));
    VidObj.FrameRate = 30;
    figure,
    set(gca, 'Color', 'k')
    set(gcf, 'Color', 'k')
    open(VidObj);
else
    figure,
end
    
for b = 1:size(SKL, 3)
    v1 = video(:, :, b)./ max(max(video(:, :, b)));
    v2 = SKL(:, :, b);
    E = bwmorph(v2, 'endpoints');
    [y1, x1] = find(E);
    dsb = sqrt((y1 - y(1)).^2 + (x1 - x(1)).^2);
    y1b = y1(dsb == min(dsb));
    x1b = x1(dsb == min(dsb));
    Db = bwdistgeodesic(SKL(:, :, b), x1b, y1b);
    v3 = v2.*0;
    v4 = v2.*0;
    v3(Db > 0 & Db < vl) = 1;
    v4(Db > ml-3 & Db < ml + 3) = 1;
    [x4, y4] = find(v4 == 1);
    x4 = mean(x4);
    y4 = mean(y4);
    frm(:,:,1) = v1+v2;
    frm(:,:,2) = v1;
    frm(:,:,3) = v1+v3;
    hold off
    imshow(frm);
    hold on
    plot(y4, x4, 'y*')
    drawnow
    pause(0.02)
    if save_stuff == 1
        this_frm = im2double(frm);
        this_frm(this_frm > 1) = 1;
        writeVideo(VidObj, this_frm);
    end
    FRM(:,:,:,b) = frm;
end
if save_stuff == 1
    close(VidObj);
    clear VidObj
end
%% segmentation of arm for data analysis
% I took out all plotting in this section

close all
nseg = 100; % divide into this number of segments
vlm = vl; % set maximum length as point selected earlier

edb = [0.2:1./nseg:1];
clear SEG

for b = 1:size(SKL,3)
    v2 = SKL(:,:,b);
    E = bwmorph(v2, 'endpoints');
    [y1, x1] = find(E);
    dsb = sqrt((y1-y(1)).^2+(x1-x(1)).^2);
    y1b = y1(dsb==min(dsb));
    x1b = x1(dsb==min(dsb));
    Db = bwdistgeodesic(SKL(:,:,b), x1b, y1b);
    Db1 = Db.*0;
    Db1(Db<=vlm) = Db(Db<=vlm);
    Dbm = Db1./max(max(Db1));
    Dbml = ml./max(max(Db1));
    Dbmn = min(find(edb > Dbml));
    clear seg
    for b1 = 1:length(edb) - 1
        clear xb yb
        Dbb = Dbm.*0;
        Dbb(Dbm > edb(b1) & Dbm < edb(b1+1)) = 1;
        [xb, yb] = find(Dbb == 1);
        seg(b1, 1) = mean(xb);
        seg(b1, 2) = mean(yb);
    end
    SEG(:,:,b) = seg;
end

%% perform calculations on segmented arm data -- segment angles
% in order to eliminate in-plane jitter along arm
% currently the metric being used for movement is the derivative of the 
% deviation that each arm segment makes from a straight line, reported 
% as the speed for that segment

clear RHO TH sym1 s1
DSS = squeeze(sqrt(diff(squeeze(SEG(:, 1, :)), 2, 2).^2 ... 
    + diff(squeeze(SEG(:, 2, :)), 2, 2).^2));
DSS = DSS - mean(DSS);

for t = 1:size(SEG,3)
    sg = squeeze(SEG(:,:,t));
    dta = squeeze(sqrt(diff(SEG(:, 1, t)).^2 ... 
        + diff(SEG(:, 2, t)).^2));
    dta1 = diff(sg);
    f1 = fspecial('average',[10 1]);
    s1(:,1) = imfilter(sg(:,1), f1,'replicate' );
    s1(:,2) = imfilter(sg(:,2), f1,'replicate' );
    f2 = fspecial('average',[5 1]);
    s2(:,1) = imfilter(sg(:,1), f2,'replicate' );
    s2(:,2) = imfilter(sg(:,2), f2,'replicate' );
    dta1 = s2 - s1;
    [theta, rho] = cart2pol(dta1(:,1), dta1(:,2));
    dd1 = length(find(dta > 20));
    dsa = abs(atan2(sin(theta), cos(theta)));
    RHO(:,t) = rho;
    TH(:,t) = dsa;  
end

% plot & save line chart
if save_stuff == 1
    figure, plot(s1(:, 2), s2(:, 1))
    hold on 
    plot(sg(:, 2), sg(:, 1))
    saveas(gcf, sprintf('%sarm_movement_line_segment_angles.png', destdirname))
    % plot & save correlation map 1
    f1 = fspecial('average',[1 5]);
    TH = imfilter(TH, f1, 'replicate' );
    %TH = TH(2:end-2,2:end-2);
    RHO = imfilter(RHO, f1, 'replicate');
    TH = abs(diff(TH,[],2));
    RHO = abs(diff(RHO,[],2));
    figure, imagesc(RHO);  %imagesc(rad2deg(abs(TH)));%rad2deg(TH));%TH(:,3:end-3))
    hold on
    plot([0 size(TH,2)], [Dbmn Dbmn],'k--','linewidth',5)
    axis square
    saveas(gcf, sprintf('%sarm_movement_correlation_map1.png', destdirname))
else
    figure, plot(s1(:, 2), s2(:, 1))
    hold on 
    plot(sg(:, 2), sg(:, 1))
    f1 = fspecial('average',[1 5]);
    TH = imfilter(TH, f1, 'replicate' );
    RHO = imfilter(RHO, f1, 'replicate');
    TH = abs(diff(TH,[],2));
    RHO = abs(diff(RHO,[],2));
    figure, imagesc(RHO);  %imagesc(rad2deg(abs(TH)));%rad2deg(TH));%TH(:,3:end-3))
    hold on
    plot([0 size(TH,2)], [Dbmn Dbmn],'k--','linewidth',5)
    axis square
end

TH = TH - min(min(TH));
TH = TH./max(max(TH));
RHO = RHO - min(min(RHO));
RHO = RHO./max(max(RHO));
sym1 = max(TH(1:Dbmn, :))./max(TH(Dbmn+1:end, :));

%% create correlation videos

if save_stuff == 1
    clear m
end
figure,
for b = 1:size(TH, 2)
    hold off
    imshow(FRM(:,:,:,b))
    hold on
    for b1 = 3:size(TH,1)
        plot(SEG(b1,2,b),SEG(b1,1,b),'c.','markersize',RHO(b1,b).*20+1)
    end
    axis([160 490 60 345])
    drawnow
    pause(0.05)
    axis off
    set(gcf,'color','k')
    m(b) = getframe;
end

hc = -1:0.05:1;
clear MN
for b = 1:size(TH,1)
    for c = 1:size(TH,1)
        c1 = corrcoef(RHO(b,1:end), RHO(c,1:end));
        xc1(b,c) = c1(3);
    end
end
        
figure, imagesc(xc1)
hold on
plot([0 size(xc1,1)], [Dbmn, Dbmn],'b--','linewidth',3)
axis equal
ind = 77;
xp = xc1(:,ind);
xp(xp < 0) = 0;
hold off
if save_stuff == 1
    saveas(gcf,sprintf('%ssegmented_arm_movement_correlation_matrix.png', destdirname))
    save(sprintf('workspace_vars_%s', date))
end

disp('continue running from live script.');
   
