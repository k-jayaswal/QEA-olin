%Converts the video file of a fidget spinner to a time signal
%by computing the average pixel value in a window for each frame
%INPUTS:
%   fname: a string of the filename of the video you want to process
%           the the video is not in the same folder, fname should include
%           the absolute path to the file
%   window_bounds: a MATLAB struct that indicates the boundaries of the
%           window to use for averaging the pixel value.
%           window_bounds.left: index of left-most pixel
%           window_bounds.right: index of right-most pixel
%           window_bounds.top: index of top-most pixel
%           window_bounds.bottom: index of bottom-most pixel
%   show_image: a boolean (0 or 1) that determines whether or not the video
%           of the fidget spinner is displayed during processing
%           set show_image to 1 if you are still trying to figure out
%           the boundaries of the window to use
%           set show_image to 0 to process the video faster
%   optional argument: you can also input whether or not to just use the
%           red, green, or blue component by passing either 
%           'red', 'r', 'green', 'g', 'blue', or 'b' as the final argument
%OUTPUTS:
%   y: a list of the averaged pixel value in the window
%   Fs: framerate of the video
function [y,Fs] = video_to_signal(fname,window_bounds,show_image,varargin)

    %If we are going to display the video, create a new figure
    if show_image
        figure();
    end
    
    %create a new VideoReader object to process the data
    vidObj = VideoReader(fname);
    
    %extract the video length and framerate
    vidDuration = vidObj.Duration;
    Fs = vidObj.FrameRate;

    %print the dimensions of the video
    disp(['Video width (pixels): ',num2str(vidObj.Width)]);
    disp(['Video height (pixels): ',num2str(vidObj.Height)]);

    %extract the window boundaries from the struct
    left = window_bounds.left;
    right = window_bounds.right;
    top = window_bounds.top;
    bottom = window_bounds.bottom;

    %make sure that the window boundaries are actually legal
    left = max(left,1);
    right = min(right,vidObj.Width);
    top = max(top,1);
    bottom = min(bottom,vidObj.Height);

    %keeps track of current frame
    count = 0;
    %used to determine when to print processing update
    prev_percent_complete = -1;

    %initialize y to a bunch of zeros with the correct length
    y = zeros(1,round(vidDuration/Fs));

    %boolean used to keep track of whether or not
    %the image of the spinner has already been created
    image_obj_exists = 0;
    %stores current frame to be displayed
    image_obj = [];

    while hasFrame(vidObj)
        %increment our fame counter at the start of each loop
        count = count+1;

        %compute progress in percent
        percent_complete = floor(100*(count/Fs) / vidDuration);

        %if the progress has incremented by 1%...
        if percent_complete ~= prev_percent_complete
            %then print the current progress
            disp(['Video processing complete: ',num2str(percent_complete),'%']);
            prev_percent_complete = percent_complete;
        end

        %extract out the next frame from the video
        vidFrame = readFrame(vidObj);
        %extract the pixel values in the desired window
        cropped_frame = vidFrame(top:bottom,left:right,:);
        
        %extract the r, g, and b values of the pixels in the desired window
        cropped_frame_red = cropped_frame(:,:,1);
        cropped_frame_green = cropped_frame(:,:,2);
        cropped_frame_blue = cropped_frame(:,:,3);

        %default: for each pixel, compute the norm of the rgb value
        v_out = uint8(sqrt(double(cropped_frame_red).^2+double(cropped_frame_green).^2+double(cropped_frame_blue).^2)/sqrt(3));

        %option 1: only use the red value
        if nargin>=7 && (strcmpi(varargin{1},'red') || strcmpi(varargin{1},'r'))
            v_out = cropped_frame_red;
        end
        
        %option 2: only use the green value
        if nargin>=7 && (strcmpi(varargin{1},'green') || lower(varargin{1},'g'))
            v_out = cropped_frame_green;
        end

        %option 3: only use the blue value
        if nargin>=7 && (strcmpi(varargin{1},'blue') || strcmpi(varargin{1},'b'))
            v_out = cropped_frame_blue;
        end

        %iterate through r, g, and b
        for n = 1:3
            %set the region in the window to be grayscale using
            %the the normalized (or only r/g/b) values
            vidFrame(top:bottom,left:right,n)=v_out;
            %visualize the window as a set of black lines in the frame
            vidFrame(max(top-1,1):min(top+1,vidObj.Height),:,n)=0;
            vidFrame(max(1,bottom-1):min(bottom+1,vidObj.Height),:,n)=0;
            vidFrame(:,max(1,left-1):min(vidObj.Width,left+1),n)=0;
            vidFrame(:,max(1,right-1):min(vidObj.Width,right+1),n)=0;
        end

        %if we are displaying the image
        if show_image
            %if we have never created the plot before
            if ~image_obj_exists
                %create a new image object
                image_obj = imshow(vidFrame);
                %set boolean to 1 (indicating the image object exists now)
                image_obj_exists = 1;

                %set tick marks
                hold on;
                axis on;
                xticks(0:100:vidObj.Width);
                yticks(0:100:vidObj.Height);

            %if we have already created the image once
            else
                %update the image object
                set(image_obj,'CData',vidFrame);
            end
            
            %show the image
            drawnow;
        end

        %compute y as the average pixel value in the window
        y(count) = sum(sum(v_out,1),2)/(size(v_out,1)*size(v_out,2));
    end
end