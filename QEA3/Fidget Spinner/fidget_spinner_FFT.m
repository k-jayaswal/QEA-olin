%uses FFT and some filtering tricks to extract angular velocity of fidget
%spinner as a function of time
%INPUTS:
%   y: list of averaged pixel values from video
%   Fs: frame rate
%   T_window: width of FFT window (in seconds) .7 is usually good
%   q: a parameter used to filter which portions of frequency curve to use
%       (filter is by width of flat portions of curve)
%       value should be in interval [0,1]
%       values closer to 1 are most "exclusive". 
%       .6 is a pretty good value to choose
%   showAnalysis: boolean that turns visualization on/off
%OUTPUTS:
%   tlist: list of times for measured values of angular frequency
%   omega_list: measured value of angular frequency
function [tlist,omega_list] = fidget_spinner_FFT(y,Fs,T_window,q,showAnalysis)
    %compute number of data points to use in FFT
    N_window = ceil(Fs*T_window);
    %perform windowed FFT of data and pull out peak frequencies
    [tlist1,freq_list1] = extract_max_integer_freq(y,Fs,N_window,showAnalysis);
    %extract the midpoints and heights of the large flat portions of signal
    [midpoint_list,height_list,width_list] = compute_contiguous_centroids(tlist1,freq_list1);
    [tlist,heights_filtererd] = filter_by_widths(midpoint_list,height_list,width_list,q);
    %de-alias signal using the fact that frequency is decreasing
    omega_list = de_alias_signal(heights_filtererd,N_window);

    %convert from integer frequencies to angular velocity
    omega_factor = 2*pi*Fs/N_window;
    omega_list = omega_factor*omega_list;

    %if showAnalysis is true, generate plots showing different stages
    %of the filtered data
    if showAnalysis
        figure();
        hold on
        plot(tlist1,omega_factor*freq_list1,'k','linewidth',1)
        plot(tlist,omega_factor*heights_filtererd,'bo','markerfacecolor','b','markersize',3);
        plot(tlist,omega_list,'ro','markerfacecolor','r','markersize',3);
        xlabel('time (sec)');
        ylabel('frequency (rad/sec)');
        title('Filtering Process');
        legend('Initial Data','Plateau Centroids','Alias Correction');
    end

    tlist = tlist-tlist(1);
end

%perform windowed FFT on data and exract the maximum
%integer frequency of the FFT for each window
%INPUTS:
%   y: list of averaged pixel values from video
%   Fs: frame rate
%   N_window: number of data points to use in each windowed FFT
%   showFFT: boolean that turns visualization on/off
%OUTPUTS:
%   tlist: list of times for each FFT
%   freq_list: maximum integer frequency of each FFT
function [tlist,freq_list] = extract_max_integer_freq(y,Fs,N_window,showFFT) 
    %used for plotting signal being processed
    tlist_base = (0:(length(y)-1))/Fs;

    %shift values of y so it has a mean of 0
    y = y-mean(y);
    %compute how many points there are in y
    Nsamples = length(y);
    %given size of window, compute the possible frequency range
    %of the FFT (using that window size)
    max_integer_freq = floor(N_window/2);
    max_integer_freq_index = max_integer_freq+1;

    %initialize the return list
    freq_list = zeros(1,Nsamples-N_window+1);
    
    %initialize plotting objects as empty arrays
    fft_fig = []; omega_fig = [];
    window_left_plot = []; window_right_plot = [];
    omega_plot = []; peak_plot = []; fft_plot = [];

    %if we are plotting, initialize the plots
    if showFFT
        %create figure showing signal and fft of window
        fft_fig=figure();
        subplot(2,1,1);
        hold on
        title('Video Signal');
        xlabel('time (sec)');
        ylabel('Average Video Color Magnitude (-)');
        plot(tlist_base,y,'k','linewidth',1);
        axis([min(tlist_base),max(tlist_base),min(y),max(y)]);
        window_left_plot = plot(0,0,'r','linewidth',1);
        window_right_plot = plot(0,0,'r','linewidth',1);
        subplot(2,1,2);
        hold on
        title('Windowed FFT of Video Signal');
        xlabel('Integer Frequency (-)');
        ylabel('Normalized FFT Magnitude (-)');
        axis([0,max_integer_freq,0,1]);
        fft_plot = plot(0,0,'k','linewidth',1);
        peak_plot = plot(0,0,'ro','markerfacecolor','r','markersize',4);

        %create figure showing peak frequencies as a function of time
        omega_fig = figure();
        hold on
        
        omega_plot = plot(0,0,'k','linewidth',1);
        axis([0,tlist_base(end),0,max_integer_freq]);
        
        title('Peak Frequency of Windowed FFT');
        xlabel('time (sec)');
        ylabel('Peak Integer Frequency (-)');
    end

    %iterate through data
    for n = 0:(Nsamples-N_window)
        %left and right indices of the window
        index0 = n+1; index1 = n+N_window;
        %data points in window
        y_vals = y(index0:index1);

        %compute and normalize fft data (divide by peak height)
        dfty = abs(fft(y_vals));
        dfty = dfty(1:max_integer_freq_index);
        dfty = dfty/max(dfty);

        %find integer frequency of peak
        [~,max_index] = max(dfty);

        %store integer frequency in output list
        freq_list(index0)=max_index-1;

        %if we are displaying the plots, update the plots
        if showFFT && mod(index0,2)==0
            %update plot showing window and corresponding FFT
            set(0,'currentfigure',fft_fig);
            set(window_left_plot,'xdata',tlist_base(index0)*[1,1],'ydata',[min(y),max(y)]);
            set(window_right_plot,'xdata',tlist_base(index1)*[1,1],'ydata',[min(y),max(y)]);
            set(fft_plot,'xdata',0:length(dfty)-1,'ydata',dfty);
            set(peak_plot,'xdata',max_index-1,'ydata',1);
            drawnow;

            %update plot showing peak frequency as function of time
            set(0,'currentfigure',omega_fig);
            set(omega_plot,'xdata',tlist_base(1:index0),'ydata',freq_list(1:index0));
            drawnow;
        end
    end
    %create tlist based on frame rate and length of freq_list
    tlist = (0:(length(freq_list)-1))/Fs;
end

%given a data set (x_i,y_i), computes the locations and width of
%the flat portions of the line plot of (x,y)
%INPUTS:
%   x: list of x coordinates of data set (assumed to be monotonic)
%   y: list of y coordinates of data set
%OUTPUTS:
%   midpoint_list: midpoint of each flat portion of data set
%   height_list: height of each flat portion of data set
%   width_list: width of each flat portion of data set
function [midpoint_list,height_list,width_list] = compute_contiguous_centroids(x,y)
    %initialize output lists
    midpoint_list = []; height_list = []; width_list = [];

    %initialize current index
    count1 = 1;
    
    %iterate until we reach end of list
    while count1<=length(y)
        %count2 is 1 + value of rightmost index of a contiguous flat region
        %initialize count2
        count2 = count1;

        %increment count2 until we reach end of flat region
        while count2<=length(y) && y(count2)==y(count1)
            count2 = count2+1;
        end
        
        %compute midpoint location and store it
        midpoint_list(end+1) = (x(count1)+x(count2-1))/2;

        %compute height and store it
        height_list(end+1) = y(count1);

        %compute width and store it
        width_list(end+1) = abs(x(count2-1)-x(count1));

        %set index to be leftmost index of next flat region
        count1 = count2;
    end
end

%given the results of compute_contiguous_centroids,
%rejects any flat candidates that are to oscillatory given their
%location in the data (as time progresses, flat portions should get wider)
%INPUTS:
%   midpoint_list: midpoints of each flat portion of data set
%   height_list: heights of each flat portion of data set
%   width_list: width of each flat portion of data set
%   q: a parameter used to filter out flat portions that are too narrow
%       value should be in interval [0,1]
%       values closer to 1 are most "exclusive". 
%       .6 is a pretty good value to choose
%OUTPUTS:
%   midpoints_filter: list of midpoints that weren't rejected
%   heights_filtered: list of heights that weren't rejected
function [midpoints_filtered,heights_filtered] = filter_by_widths(midpoint_list,height_list,width_list,q)
    %list of max width value seen up to this point
    max_width_list = 0*width_list;
    %initialize max_width_list
    max_width_list(1) = width_list(1);

    %populate max_width_list
    for n = 2:length(width_list)
        max_width_list(n)=max(max_width_list(n-1),width_list(n));
    end

    %initialize return lists
    midpoints_filtered = []; heights_filtered = [];

    %iterate through data and reject and flat regions in data that are
    %not sufficiently wide given placement in data.
    for n = 1:length(width_list)
        if width_list(n)>=q*max_width_list(min(max(n,4),length(max_width_list)))
            midpoints_filtered(end+1)=midpoint_list(n);
            heights_filtered(end+1)=height_list(n);
        end
    end
end

%uses the fact that the angular frequency of the fidget spinner is
%monotonically decreasing to back out the correct peak frequency
%even, when the fidget spinner is rotating multiple times per frame
%INPUTS:
%   y: list of peak integer frequencies from windowed FFT
%       (already filtered by looking at flat portions)
%   N_fft_samples: number of sample points used in FFT
%OUTPUTS:
%   y: "adjusted" list of peak integer frequencies
function y = de_alias_signal(y,N_fft_samples)
    %iterate backwards through each element of y
    for n = (length(y)-1):-1:1
        %use aliasing formula to compute the lowest integer frequency
        %that would alias to y(n), but would still have a larger value
        %that y(n+1) (assuming that y(n+1) has already been corrected

        temp1 = y(n);
        temp1 = temp1+ N_fft_samples*ceil((y(n+1)-temp1)/N_fft_samples);
        while temp1<=y(n+1)
            temp1 = temp1+N_fft_samples;
        end

        temp2 = -y(n);
        temp2 = temp2+ N_fft_samples*ceil((y(n+1)-temp2)/N_fft_samples);
        while temp2<=y(n+1)
            temp2 = temp2+N_fft_samples;
        end

        y(n) = min(temp1,temp2);
    end
end