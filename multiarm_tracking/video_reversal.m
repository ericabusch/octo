%% multiple arm tracking script
% script for tracking multiple arm movement in a video frame
% takes as input a video featuring a whole octopus body
% select a key frame where all arms are visible. Work forward and backward
% from this key frame to identify where arms are crossed.
% ELB

clear all
save_stuff = 0; % do you want to save the new video generated?
get_input = 1;
key_frame_time = 3.0; % at what second in the video is your key frame? 
times_cross = [13.0,14.0]; 

if get_input == 1
    clear key_frame_time times_cross
    
    prompt1 = 'Input time of key frame, when all arms are visible: ';
    key_frame_time = str2double(input(prompt1, 's'));
    
    prompt2 = 'Input time when arms cross or become obscured within []: ';
    times_cross = input(prompt2);

end

disp("Select your original raw .avi file");
[raw_vid_name, raw_vid_loc] = uigetfile('*.avi');
cd(raw_vid_loc);

if save_stuff == 1
    save_path = sprintf('fwd_bkwd_%s_', raw_vid_name);
end

sbs1 = 2; % spatial subsampling

% read in the orig video
v = VideoReader(raw_vid_name);
% v = VideoReader(sprintf('%s_%s',raw_vid_loc, raw_vid_name));
fps = v.FrameRate;
p = 0;
frames_crossed = [];

% go through each frame
while hasFrame(v)
    % if the current frame is the identified key frame, set it as such
    if v.CurrentTime == key_frame_time
        kf = readFrame(v);
    end
    
    % if thee current frame is a frame where the arms are crossed,
    % add this frame to a list tracking those frames
    if ismember(v.CurrentTime, times_cross)
        frames_crossed(end+1) = p; %#ok<*SAGROW>
    end
        
    p = p+1;
    video1 = readFrame(v);
    video(:,:,p) = double(video1(1:sbs1:end, 1:sbs1:end, 1));
end


for f = 1:size(frames_crossed,2)
    fprintf('size frames_crossed is %d; f is %d\n',size(frames_crossed), f);
    clear new_vid;
    
    new_vid = VideoWriter(sprintf('%s_%d', raw_vid_name, f)); %#ok<TNMLP>
    new_vid.FrameRate = 30;
    fprintf('%s_%d\n', raw_vid_name, f)
    x = frames_crossed(f);
    fprintf('running frame: %d\n',x)

    open(new_vid);
    
    for y = 1:x+10
        if y <= x            
            writeVideo(new_vid, uint8(video(:,:,y)));
        else
            z = x - (y - x);
            writeVideo(new_vid, uint8(video(:, :, z)));
        end
    end
    fprintf('done frame: %d\n',x)

    close(new_vid);
    fprintf('closed',x);

end








