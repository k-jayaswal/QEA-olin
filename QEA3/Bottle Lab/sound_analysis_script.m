function sound_analysis_script()
    % Collects sound data and then analyzes it.
    % Usage:
    %   Run the script; it will enter recording mode for 5 seconds.
    %   While in recording mode, play the sound to analyze.
    %     Ensure the sound source is close enough to your computer to be
    %     heard.
    %   After recording has ended, the script will determine the frequency
    %   that has the *largest* DFT coefficient. This is the Primary
    %   Frequency Component, which will be printed to your Command Window.
    %     To help you determine whether this frequency component matches
    %     the sound you played, the script will also play back a pure tone
    %     at the Primary Frequency Component. You can tell by ear if the
    %     script picked up what you intended (i.e., you can tell by ear if
    %     the Primary Frequency Component is very different from the sound
    %     you just recorded.)

    run_sound_test();

    %extra script I wrote to make sure I understood how FFT works
    % fft_test();
end

%this function does the following
%1. collects a sound sample
%2. computes the highest amplitude frequency
%      (within a specified frequency band)
%3. plays back the primary frequency component for validation
function run_sound_test()
    %set the sampling rate (in Hz)
    % Fs = 8000;
    % Fs = 48000;
    Fs = 192000;

    %set the number of bits per sample
    % BitsPerSample = 8;
    % BitsPerSample = 16;
    BitsPerSample = 24;

    %set number of channels (mono or stereo)
    % NumChannels = 2;
    NumChannels = 1;

    %create an audiorecorder object
    recObj = audiorecorder(Fs, BitsPerSample, NumChannels);

    %set the amount of time to record (in seconds)
    recDuration = 5;

    disp(['Recording Has Begun, Duration: ',num2str(recDuration),' Seconds']);

    %record sound from soundcard (result is stored in recObj)
    recordblocking(recObj,recDuration);

    disp("Recording Has Ended");

    %uncomment if you'd like to play back the recorded sound
    % play(recObj);

    %store sound signal in a vector
    y = getaudiodata(recObj);

    %generate vector of points for each sample time (in seconds)
    t_list = (0:length(y)-1)/Fs;

    min_freq = 300; %min cutoff frequency in rad/sec 
    max_freq = 1e5; %max cutoff frequency in rad/sec

    %compute primary frequency component of recorded sound
    %that is within the desired frequency band
    omega_primary = estimate_freq_from_fft(t_list,y,min_freq,max_freq);

    %print out frequency of sound in rad/sec
    disp(['Primary Frequency Component: ',num2str(omega_primary),' rad/sec']);

    %generate a sine wave at the primary frequency
    pure_sound = 3*sin(omega_primary*t_list);

    %uncomment to play back pure tone
    sound(pure_sound,Fs);
end

%computes primary frequency component of signal
%(within a specified frequency band)
%INPUTS:
%t_list: list of sampling time-points (assumed in seconds)
%y: measured signal (units irrelevant)
%min_cutoff_freq: minimum cutoff frequency (any frequency component below
%                   this threshold is ignored)
%max_cutoff_freq: maximum cutoff frequency (any frequency component above
%                   this theshold is ignored)
function omega_out = estimate_freq_from_fft(t_list,y,min_cutoff_freq,max_cutoff_freq)
    %compute fft of signal
    dft_signal = abs(fft(y));

    %compute number of sample points
    num_sample_points = length(y);

    %compute the sampling frequency (in Hz). Note that in the interval between
    %t(1) and t(end), num_sample_points-1 periods have occured
    %(number of fences is equal to number of posts-1)
    %hence the -1 in the numerator
    Fs = (num_sample_points-1)/(t_list(end)-t_list(1));

    %compute list of integer frequencies
    omega_list_integer = 0:(num_sample_points-1);

    %ratio of integer frequency to rad/sec frequency
    int_to_radsec_ratio = (2*pi*Fs/num_sample_points);

    %generate list of frequencies (in rad/sec)
    omega_list = omega_list_integer*int_to_radsec_ratio;

    % Convert frequencies to Hz
    % NOTE: You can plot using this as your x-axis (in place of omega_list)
    % to report frequencies in Hz (rather than rps)
    freq_list = omega_list / (2 * pi);

    %compute the maximum integer frequency of fourier basis used in the FFT
    max_integer_freq = floor(num_sample_points/2);

    %compute the corresponding maximum index of 
    % fourier basis used in the FFT
    max_index = max_integer_freq+1;

    min_index = 1;

    %compute the integer frequencies for the min/max cutoff freq
    min_cutoff_integer_freq = floor(min_cutoff_freq/int_to_radsec_ratio);
    max_cutoff_integer_freq = ceil(max_cutoff_freq/int_to_radsec_ratio);

    %compute dft index for min/max cutoff_freq
    min_cutoff_index = min_cutoff_integer_freq+1;
    max_cutoff_index = max_cutoff_integer_freq+1;

    %keep indices within correct index range of 1 to floor(num_sample_points/2)
    index_left = min(max(min_index,min_cutoff_index),max_index);
    index_right = min(max(min_index,max_cutoff_index),max_index);

    %compute the amplitude and index of the primary frequency component
    %(the primary frequency has the maximum amplitude in desired freq band)
    [output_amplitude,temp_index]= max(dft_signal(index_left:index_right));
    output_index = temp_index+(index_left-1);

    %compute primary integer frequency using the index
    output_integer_freq = output_index-1;

    %compute the primary frequency in rad/sec
    omega_out = output_integer_freq*int_to_radsec_ratio;

    %for validation, generate a sinusoid function that will give fft with
    %a single spike located at the primary frequency (in the band)
    %and have the sample amplitude as in the original signal
    test_sinusoid = @(t_in) 2*(output_amplitude/num_sample_points)*cos(omega_out*t_in);

    %evaluate sinusoid across the sampling time points
    y_test = test_sinusoid(t_list);

    %compute the dft of the validation sinusoid
    dfty_test = abs(fft(y_test));

    %limit the plotting window to frequency of relevant amplitudes
    max_plot_index = max_index;
    while dft_signal(max_plot_index)<output_amplitude/100
        max_plot_index = max_plot_index-1;
    end

    %plot fft of signal, and compare with fft of validation sinusoid.
    figure();
    hold on;
    plot(omega_list(1:max_plot_index),dft_signal(1:max_plot_index)/(num_sample_points),'b');
    plot(omega_list(output_index),dfty_test(output_index)/(num_sample_points),'ro','markerfacecolor','r','markersize',5);

    plot(omega_list(1:(index_left-1)),dft_signal(1:(index_left-1))/(num_sample_points),'k');
    plot(omega_list((index_right+1):max_plot_index),dft_signal((index_right+1):max_plot_index)/(num_sample_points),'k');
    
    xlabel('frequency (rad/sec)');
    ylabel('abs(FFT(signal)) (units unclear)');
    title('Fourier analysis of recorded sound')
    legend('Desired Band','Frequency selected')
end

%This is a check to make sure that I properly understand the indexing
%and overall functionality of MATLAB's FFT
%To do this, I generate a phase shifted sine wave with a frequency that
%is specifically chosen to generate a single spike for the FFT
%then, I generate the predicted FFT and compare the two
function fft_test()
    
    %choose a random number of sampling points
    num_sample_points = randi(30)+5;

    %choose a random sampling frequency
    Fs = exp(3*randn());

    %generate list of sampling times
    t_list = (0:(num_sample_points-1))/Fs;

    %randomly choose the integer frequency of sine wave
    %bias towards either constant frequency or
    %maximum frequency of fourier basis used in this FFT
    choice = randi(3);
    q = 0;
    if choice==1
        disp('choosing constant function');
        q = 0;
    end
    if choice==2
        disp('choosing maximum frequency of fourier basis');
        q = floor(num_sample_points/2);
    end
    if choice==3
        disp('choosing random frequency (except zero)');
        q = randi(floor(num_sample_points/2));
    end

    %compute frequency in radians/time
    omega = (2*pi*Fs)*(q/num_sample_points);

    %choose a random phase shift and amplitude
    phi = 2*pi*rand();
    A = exp(2*randn());

    %generate the sine wave function
    sine_wave_func = @(t_in) A*cos(omega*t_in+phi);

    %evaluate sine wave function at time points
    y = sine_wave_func(t_list);

    %compute fft of signal
    dfty = fft(y);

    %generate list of frequencies (in rad/sec)
    omega_list = (0:(num_sample_points-1))*(2*pi*Fs/num_sample_points);

    %compute predicted FFT
    dfty_predicted_real = zeros(length(y),1);
    dfty_predicted_imag = zeros(length(y),1);
    
    
    index1 = q+1;
    index2 = length(y)+1-q;
    if index2>length(y)
        index2 = index2-length(y);
    end

    if q==0 || q==num_sample_points/2
        %if frequency was zero, or the integer frequency was exactly half 
        %of the number of sample points, then there is no imaginary
        %component in this case, and amplitude is doubled
        dfty_predicted_real(index1) = A*cos(phi)*length(y);
    else
        %compute predicted dft in regular case
        %two spikes with cojugate symmetry
        dfty_predicted_real(index1)=A*cos(phi)*length(y)/2;
        dfty_predicted_imag(index1)=A*sin(phi)*length(y)/2;
        dfty_predicted_real(index2)=dfty_predicted_real(index1);
        dfty_predicted_imag(index2)=-dfty_predicted_imag(index1);
    end

    %plot comparison of actual and predicted FFT
    hold on
    plot(omega_list,real(dfty),'k','linewidth',2);
    plot(omega_list,imag(dfty),'b','linewidth',2);
    plot(omega_list,dfty_predicted_real,'r--','linewidth',2);
    plot(omega_list,dfty_predicted_imag,'g--','linewidth',2);
    xlabel('frequency (rad/time)');
    ylabel('FFT(y) (-)');
    title('Testing FFT on a sine wave with frequency chosen to make a spike');
    legend('re(FFT)','im(FFT)','re(Predicted FFT)','im(Predicted FFT)');
end