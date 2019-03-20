%function that overlaps the image with visual representation 

function image_overlap = image_overlapping(seg_coordinates, original_data, mapped_parameters, model_type)
    %Find the axis of where the layer of segmentation is
    
    num_coords = length(seg_coordinates); 
    
    x_diff = seg_coordinates(1, 1) - seg_coordinates(1, num_coords); 
    y_diff = seg_coordinates(2, 1) - seg_coordinates(2, num_coords);
    z_diff = seg_coordinates(3, 1) - seg_coordinates(3, num_coords); 
    
    %segmentation axis is 1 = x-axis, 2 = y-axis and 3 = z-axis 
    seg_axis = 1; 
    
    if (y_diff == 0)
        seg_axis = 2; 
    elseif (z_diff == 0)
        seg_axis = 3; 
    end 
    
    %now we get the anatomical image from the original data
    %This is the configuration to check the correct orientation of the
    %image
    if (seg_axis == 1) %this first one is corrected 
        canvas = squeeze(original_data(seg_coordinates(1,1), :, :));
    elseif (seg_axis == 2) 
        canvas = squeeze(original_data(seg_coordinates(2,1), :, :)); 
    else
        canvas = squeeze(original_data(:, : ,seg_coordinates(3,1))); 
    end 
    
    %next we check if we can hide the segmentation from the canvas 
    
    %there shoould be an if statement checking for which axis we are on,
    %but for now we just know we are on the x-axis
    
    for coordinate = 1 : num_coords
        %take the x coordinate as the z 
        x_val = seg_coordinates(3, coordinate); 
        y_val = seg_coordinates(2, coordinate); 
        canvas(y_val, x_val) = 0; 
    end 
    
 
    canvas = imadjust(uint8(mat2gray(canvas)*255));
    
    
    %Here we switch the type of model that we want 
    %1 = adc mapping 
    %2 = best model mapping
    %3 = directional mapping 
    
    %set up a canvas 1 for each model
    [canvas_rows, canvas_cols] = size(canvas); 
    
    lots_of_canvas = zeros(canvas_rows, canvas_cols, 3, length(mapped_parameters(1, :, 1)));
    for models = 1 : 13
        for RGB_layer = 1 : 3
            lots_of_canvas(:, :, RGB_layer, models) = canvas; 
        end 
    end
        
    
    switch (model_type)
        case 1 
            %user has selected ADC mapping! 
            %loop through every different model from the paramters
            for model_no = 1 : 13
                
                %loop through all the voxels we have in mapped_parameters 
                for voxel_no = 1 : length(mapped_parameters(:, 1, 1));
                    %we are taking the z value to be the x value here
                    x_val = mapped_parameters(voxel_no, model_no, 16); 
                    y_val = mapped_parameters(voxel_no, model_no, 15); 
                    adc_val = double(mapped_parameters(voxel_no, model_no, 1)); 
                    max_adc_val = double(max(mapped_parameters(:, model_no, 1)));
                    %Now we need to produce the colour gradient for the adc
                    
                    max_red = 66; 
                    min_red = 244;
                    
                    %first we normalise the voxel's adc with the maximum
                    normalised_adc = adc_val / max_adc_val; 
                    green_val = uint8(round(min_red - (min_red*normalised_adc), 0));
                    if(green_val < max_red)
                        green_val = max_red;
                    end
                    %populate canvas with the colour on three dimensions 
                    red_val = 244;
                    blue_val = 66;
                    lots_of_canvas(y_val, x_val, 1, model_no) = red_val; 
                    lots_of_canvas(y_val, x_val, 2, model_no) = green_val; 
                    lots_of_canvas(y_val, x_val, 3, model_no) = blue_val;
                end
            end 
   
            
        case 2
            %user has selected best model mapping! 

            %loop through each voxel
            for voxel_no = 1 : length(mapped_parameters(:, 1, 1))
                %we are taking the z value to be the x value here
                x_val = mapped_parameters(voxel_no, 1, 16); 
                y_val = mapped_parameters(voxel_no, 1, 15); 
                %get the best fitted model by getting the index of the minimum BIC 
                min_BIC = min(mapped_parameters(voxel_no, :, 13)); 
                min_BIC_index = find(mapped_parameters(voxel_no, :, 13) == min_BIC);
                
                %set colour code for all 13 models 
                red_val = 0;
                green_val = 0;
                blue_val = 0; 
                
                switch (min_BIC_index)
                    case 1 %red
                        red_val = 244;
                        green_val = 66;
                        blue_val = 66; 
                    case 2 %orange 
                        red_val = 244;
                        green_val = 164;
                        blue_val = 66;
                    case 3 %yellow 
                        red_val = 232;
                        green_val = 244; 
                        blue_val = 66;
                    case 4 %brown 
                        red_val = 140; 
                        green_val = 72; 
                        blue_val = 5; 
                    case 5 %Green 
                        red_val = 167; 
                        green_val = 244; 
                        blue_val = 66;
                    case 6 %Dark.G
                        red_val = 27; 
                        green_val = 109; 
                        blue_val = 37; 
                    case 7 %cyan 
                        red_val = 66; 
                        green_val = 244; 
                        blue_val = 209; 
                    case 8 %blue 
                        red_val = 66; 
                        green_val = 173; 
                        blue_val = 244; 
                    case 9 %navy
                        red_val = 35;
                        green_val = 37; 
                        blue_val = 109;
                    case 10 %purple 
                        red_val = 203; 
                        green_val = 66;
                        blue_val = 244;
                    case 11 %pink 
                        red_val = 244; 
                        green_val = 66;
                        blue_val = 200; 
                    case 12 %lightpink 
                        red_val = 249; 
                        green_val = 199;
                        blue_val = 245; 
                    case 13 %mauv 
                        red_val = 96;
                        green_val = 10; 
                        blue_val = 89;
                end
               %colour has been set, now we set the colour  
                lots_of_canvas(y_val, x_val, 1, 1) = red_val; 
                lots_of_canvas(y_val, x_val, 2, 1) = green_val;
                lots_of_canvas(y_val, x_val, 3, 1) = blue_val; 

            end    
    end
    
    
    

    switch model_type
        case 1 
            %display all the models mappings
            figure;
            for model_no = 1 : 13
                subplot(4, 4, model_no);
                imshow(uint8(squeeze(lots_of_canvas(:,:,:,model_no)))); 
                title(model_no);
                image_overlap = lots_of_canvas; 
            end 
            
        case 2 
            figure;
            subplot(1,1,1);
            image_overlap = uint8(squeeze(lots_of_canvas(:,:,:,1)));
      	    imshow(image_overlap); 
            title("RAINBOW");
           
    end
end 