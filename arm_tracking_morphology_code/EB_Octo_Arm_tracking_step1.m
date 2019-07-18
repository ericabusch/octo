%% sequential video processing 
% eb copy


%% Import video and background image
% notes for input: the task video should be included as .avi and a single 
% image of background as a .png. 
% Both of these are necessary to start the program. 

imp = 1; %import new file?
% open the video
if imp == 1
    clear all
end
make_changes = 0; % are we going to have to clip the video?
save_stuff = 1; % do we want to save the output?
if save_stuff == 1
    if make_changes == 1
        destdir = sprintf('NEWskeletonized_arm_videos_%s', date);
    else
        destdir = sprintf('skeletonized_arm_videos_%s', date);
    end
    
    if exist(destdir, 'dir') == 7
        rmdir(destdir)
    end
    mkdir(destdir)
end
 
sbs1 = 2; %spatial sub-sampling
% import raw avi file
disp("Select raw .avi file.");
[fn, dn] = uigetfile('*');
cd(dn)

%132-
root_dir = cd;
v = VideoReader(fn);
p = 0;
while hasFrame(v)
    p = p + 1;
    video1 = readFrame(v);
    video(:,:,p) = double(video1(1:sbs1:end, 1:sbs1:end, 1));
end

% Get background image
% clip background image from first fraame, then trim video
if make_changes == 1
    BK = video(:,:,1);
    video = video(:,:,900:1200);
else
    % choose background image file and proceed.
    disp("Select raw background image file.");
    [bk,dn] = uigetfile('*.png');
    BK = double(imread(sprintf('%s/%s',dn,bk))); %read in a background image
    BK = BK(1:sbs1:end, 1:sbs1:end, 1);
end
%% define the variables
kf = 1; % define keyframe
BWTH = 56; % BW threshold
bright = 0; % is the animal bright?
if bright == 1
    BWTH = 80;
end
dsk1 = 3; % disk element diameter
GS1 = fspecial('average',[5 5]);
arm1 = 30; %arm diameter for initial tracking
arm_V1 = 50; %arm diameter during tracking (should be much larger)

%% skeletonization of the arm
figure,
clear ST
for t = 1:size(video,3)-1
    if t > 1 
        s_test = sum(sum(bwl));
    end
    image1 = video(:,:, kf + t); 
    % apply filter to the frame
    image1 = imfilter(image1, GS1);
    imageb = image1;
    if bright == 0
        image1 = abs(image1 - BK);
    else
        image1 = (image1 - BK);
    end

    if bright == 1
        BW = image1;
        BW(BW > BWTH) = 1;
        BW(BW ~= 1) = 0;
    else
        BW = image1;
        BW(BW < BWTH) = 0;
        BW(BW ~= 0) = 1;
    end
    % label objects
    BWL = bwlabel(BW);
    ip = unique(BWL(BWL ~= 0));
    BWL1 = BWL.*0;

    if length(ip) > 1
        clear ls LS
        for s = 1:length(ip)
            ls(s) = length(find(BWL == ip(s)));
        end
        LS = ip(ls == max(ls));
        LS = LS(1);
        BW(BWL ~= ip(LS)) = 0;
        se = strel('disk', dsk1);
        BWn = imdilate(BW, se); %potentially not doing anything
    end
    
    BW = imdilate(BW,se);
    BWL1(ismember(BWL, ip(ls > 25))) = 1;
    DS1 = bwdist(1 - BW); %distance to nearest black pixel
    [mdx, mdy] = find(DS1 == max(max(DS1))); %center of mass (of the mantle)
    SK = bwmorph(BW, 'skel', Inf);
    TL = SK.*DS1;
    TL1 = TL.*0;
    TL2 = TL1.*0;
    TL1(TL>0 & TL<1000) = 1;
    TL2(TL>0 & TL<arm1) = 1;
    bwl = bwlabel(TL2);
    STATS = regionprops(bwl, 'Area');
    clear areaa
    for p = 1:length(STATS)
        areaa(p) = STATS(p).Area;
    end
    arms1 = find(areaa < 25);
    bwl(ismember(bwl, arms1)) = 0;
    BT = bwl.*0;
    BT(bwl > 0) = 1;
    bwl = bwlabel(BT);
    hold off
    imagesc(bwl)
    drawnow
    pause(0.01)
    s_t = sum(sum(bwl));
    if t > 1
        ST(kf + t + 1) = s_t - s_test;
    end
    SKL(:,:,t) = bwl;
end

%% output from first step is skeletonized version of the video. 
% If this is acceptable, save it here for later analysis.

clear m
VidObj = VideoWriter(sprintf('%s/skeleton1.avi', destdir));
VidObj.FrameRate = 30;
 
 figure,
 set(gca,'Color','k')
 set(gcf,'Color','k')
open(VidObj);

for b = 1:size(SKL, 3)
    imshow(SKL(:,:,b));
    writeVideo(VidObj,SKL(:,:,b));
end

close(VidObj);
clear VidObj

save(sprintf('%s/skeleton1.mat', destdir),'SKL')
fprintf('saved in %s\n', sprintf('%s/skeleton1.mat', destdir));
