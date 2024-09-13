global frame_id;
global no_of_channels;
global int_size;
global chirp_size;
global no_of_rows;
global frame_size;
global file_name;
global numADCSamples;
global maxADCIndex;

global velocityOfLight;
global sampleRate;
global freqSlope;
global carrierFrequency;

frame_id = 1000;
no_of_channels = 4;
int_size = 16;
chirp_size = 512;    %256 complex numbers not bytes
no_of_rows = 128;
frame_size = chirp_size * no_of_rows;     % 2 for 2 integer values in a complex number = 1048576
file_name = 'D:\\Radar\\From LRDE\\Moving1i1o.bin';
numADCSamples = 256;
maxADCIndex = 128;

velocityOfLight = physconst('LightSpeed');
sampleRate = 5000e3;
freqSlope = 8e12;
carrierFrequency = 77e9;


[max_no_frames, filesize] = get_max_frames();
[result, lvds_data] = read_file();
i_and_q_data = populate_channel_data(lvds_data);
processed_i_and_q = sum_channel_data(i_and_q_data);

abs_i_and_q = populate_abs_data(processed_i_and_q, true);
plot_surface_chart(abs_i_and_q, 121, "x-axis", "y-axis", "RAW I and Q data");

one_dim_fft = perform_1D_fft(processed_i_and_q);
abs_fft = populate_abs_data(one_dim_fft, false);
plot_surface_chart(abs_fft.', 122, "x-axis", "y-axis", "One Dim FFT");

two_dim_fft = perform_2D_fft(processed_i_and_q);
abs_fft = populate_abs_data(two_dim_fft, true);
plot_surface_chart(abs_fft, 123, "Time", "Speed", "Two Dim FFT");

function [max_no_frames, filesize] = get_max_frames()
    global file_name;
    global frame_size;
    s = dir(file_name);
    filesize = s.bytes;
    max_no_frames = cast(filesize / frame_size, "uint32");
end

function [result, lvds_data] = read_file()
    global file_name;
    fid = fopen(file_name,'r');
    lvds_data = fread(fid, 'int16');
    result = true;
end

function  adcData = populate_channel_data(lvds_data)
    global numADCSamples;
    global no_of_channels;
    global numChirps;
    filesize = size(lvds_data, 1);
    numChirps = filesize/2/numADCSamples/no_of_channels;
    i_and_q_data = zeros(1, filesize/2);
    %combine real and imaginary part into complex data
    %read in file: 2I is followed by 2Q
    counter = 1;
    for i=1:4:filesize-1
        i_and_q_data(1,counter) = lvds_data(i) + sqrt(-1)*lvds_data(i+2); 
        i_and_q_data(1,counter+1) = lvds_data(i+1)+sqrt(-1)*lvds_data(i+3);
        counter = counter + 2;
    end
    % create column for each chirp
    i_and_q_data = reshape(i_and_q_data, numADCSamples*no_of_channels, numChirps);
    %each row is data from one chirp
    i_and_q_data = i_and_q_data.';

    adcData = zeros(no_of_channels,numChirps*numADCSamples);
    for row = 1:no_of_channels
     for i = 1: numChirps
        adcData(row, (i-1)*numADCSamples+1:i*numADCSamples) = i_and_q_data(i, (row-1)*numADCSamples+1:row*numADCSamples);
     end
    end
end

function processed_i_and_q = sum_channel_data(adc_data)
    global numADCSamples;
    global no_of_channels;
    global numChirps;
    global maxADCIndex;
    global frame_id;
    global frame_size;

    colNum = frame_id * numADCSamples;
    processed_i_and_q = zeros(maxADCIndex, numADCSamples);
    i_and_q_data = zeros(1, maxADCIndex*numADCSamples);
    for rowIndex = 1:maxADCIndex
        for colIndex  = 1:numADCSamples 
            i_and_q_data(colNum) = adc_data(1,colNum) + adc_data(2,colNum) + adc_data(3,colNum) + adc_data(4,colNum);
            colNum = colNum +1;        
        end
    end
    col_num = (frame_id -1)* numADCSamples +1;
    for row_index = 1:maxADCIndex
        for col_index  = 1:numADCSamples 
            processed_i_and_q(row_index,col_index) = i_and_q_data(col_num);
            col_num = col_num +1;
        end
    end
end

function one_dim_fft = perform_1D_fft(processed_i_and_q)
    global numADCSamples;
    global maxADCIndex;
    one_dim_fft = fft(processed_i_and_q); %,maxADCIndex,numADCSamples);
end
    
function two_dim_fft = perform_2D_fft(processed_i_and_q)
    global numADCSamples;
    global maxADCIndex;
    two_dim_fft = fft2(processed_i_and_q,maxADCIndex,numADCSamples);
end

function abs_fft = populate_abs_data(fft_array, isTwoDimFFT)
    global numADCSamples;
    global maxADCIndex;

    if isTwoDimFFT == true
        abs_fft = 20*log10(abs(fft_array));
    else
        abs_fft = 10*log10(abs(fft_array));
    end
end

function plot_surface_chart(data, chart_id, x_label, y_label, chart_title)
    global frame_id;
    figure(chart_id)
    surf(data);
    view(0,90)
    axis tight
    shading flat;
    colormap jet;
    view(0,90);
    xlabel(x_label);
    ylabel(y_label);
    frame_id_text = "Frame Number " + frame_id;
    title(chart_title,frame_id_text);
    colorbar;
    grid on;
end

function plot_image(rngdopresp,rnggrid,dopgrid, chartId)
    global frame_id;
    figure(chartId);
    rngdopresp = flip(rngdopresp,1);
    imagesc(rnggrid,dopgrid,20*log10(abs(rngdopresp))')
    ylabel('Relative Speed (m/s)')
    xlabel('Range (m)')
    axis xy;
    colorbar;
    colormap jet;
    frame_id_text = "Frame Number " + frame_id;
    title('Range Vs Speed Map',frame_id_text);
    grid on;
end
