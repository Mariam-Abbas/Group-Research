%Function for calculating the adc of each voxel of a segment from the data 
function adc_matrix = adc_segment(segment_file)
    global data; 
    global trueSignal;
    %Identify the coordinates of the relevant voxels from the segment 
    voxel_coordinates = mask_coordinates_extractor(segment_file); 
    
    %get the size ofthe coordinate matrix 
    [m, n] = size(voxel_coordinates);
    
    %For every voxel we calculate 
    % (i) The true signal 
    % (ii) The minimised adc and s0 value
    
    %initilise an array to store the adc and s0 values with its coordinates
    voxel_property = zeros(5, n); %where n is the number of voxel entries 

    for voxel = 1 : n 
        %find the x y and z coordinates of the voxel from the matrix 
        voxel_x = voxel_coordinates(1, 1); 
        voxel_y = voxel_coordinates(2, 1); 
        voxel_z = voxel_coordinates(3, 1); 
        %caluclate the true signal 
        trueSignal = squeeze(data(voxel_x, voxel_y, voxel_z, :));
        %ready the fmincon function );
disp(ADCopt);
disp(best_ADC);
plot_ADC(best_ADC);

%Function for calculating the adc of each voxel of a segment from the data 
function adc_matrix = adc_segment(segment_file)
    global data; 
    global trueSignal;
    %Identify the coordinates of the relevant voxels from the segment 
    voxel_coordinates = mask_coordinates_extractor(segment_file); 
    %how paddy will do this: take global data -> take seg_file -> loop through each dimension, get an if statement  
    
    %get the size ofthe coordinate matrix 
    [m, n] = size(voxel_coordinates);
    
    %For every voxel we calculate 
    % (i) The true signal 
    % (ii) The minimised adc and s0 value
    
    %initilise an array to store the adc and s0 values with its coordinates
    voxel_property = zeros(5, n); %where n is the number of voxel entries 

    for voxel = 1 : n 
        %find the x y and z coordinates of the voxel from the matrix 
        voxel_x = voxel_coordinates(1, 1); 
        voxel_y = voxel_coordinates(2, 1); 
        voxel_z = voxel_coordinates(3, 1); 
        %caluclate the true signal 
        trueSignal = squeeze(data(voxel_x, voxel_y, voxel_z, :));
        %ready the fmincon function 
        lb = [0 0];    
        x0 = [0.01; 0.01];
        [optADC, optS0] = fmincon(@cominedOptimise, x0, [], [], [], [], lb, []);
        %insert the minimised ADC and S0 into the voxel property matrix 
        voxel_property(4, voxel) = optADC; 
        voxel_property(5, voxel) = optS0; 
        %insert the coordinates of the voxel into its property too 
        voxel_property(1, voxel) = voxel_x; 
        voxel_property(2, voxel) = voxel_y; 
        voxel_property(3, voxel) = voxel_z; 
    end
    
    
    %return the adc properties of each voxel 
    adc_matrix = voxel_property; 
    
    %beyond here, the code will throw a greyscaled image of the adc of the
    %segment 
end
