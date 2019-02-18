%Extraction of the coordinates of the segmentation mask 
function [layers, coordinates] = mask_coodinates_extractor(filename)
    %Load the segmentation file 
    seg_data = double(niftiread(filename));
    [m,n,o] = size(seg_data)
    %loop through the entire segmentation the first time to find storage
    %parameters 
    
    %We need to hold the number of layers and the number of coordinates to
    %preset our array sizes when we populate them in the second loop 
    num_layers = 0;
    num_coordinates = 0; 
    for z = 1 : o
        %Keep track of when we have found a non_zero 
        non_zero_found = 0;
        for y = 1 : m
            for x = 1 : n 
                if seg_data(x, y, z) > 0
                    non_zero_found = 1; 
                    num_coordinates = num_coordinates + 1; 
                end
            end
        end
        %If a non_zero has been found, we know this layer contains part of
        %the segment
        if non_zero_found == 1
            num_layers = num_layers + 1; 
        end
    end

    %initialise the matrix storage for all the coordinates and a second matrix
    %for just storing the list of layers 
    coordinate_matrix = zeros(3, num_coordinates);
    layer_array = zeros(1, num_layers);
    %hit_num keeps track of number of hits we have got to increment our matrix
    hit_num = 0;
    layer_hits = 0;
    %loop the second time to extract the coordinates
    for z = 1 : o
        %keep track of whether we got a hit
        hit_found = 0; 
        for y = 1 : m
            for x = 1 : n 
                if seg_data(x, y, z) > 0
                    %When we find a non_zero incremenet the number of hits
                    hit_num = hit_num + 1;
                    %insert x, y and z into the coordinate matrix
                    coordinate_matrix(1, hit_num) = x;
                    coordinate_matrix(2, hit_num) = y; 
                    coordinate_matrix(3, hit_num) = z; 

                    %we found a hit, set hit_found to true
                    hit_found = 1; 
                end
            end
        end
        %If hit was found we can add this layer to the layer array
        if hit_found == 1 
            %incremenet the layer hit index
            layer_hits = layer_hits + 1;
            layer_array(1, layer_hits) = z; 
        end
    end
    %return the values 
    coordinates = coordinate_matrix;
    layers = layer_array;
end

