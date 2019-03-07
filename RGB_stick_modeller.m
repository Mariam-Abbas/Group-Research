%function that obtains the RGB values required to graphically model a stick
%model that is fitted 
function RGB_stick_model = RGB_stick_modeller(normalised_adc, theta, phi, grid_size)
    %First we need to identify the coorindates of the stick model from its 
    %angle orientation 


    %identify the length of the stick model first 
    %Assumption: theta will be the x plane and phi will be the y plane 
    y_value = grid_size*sin(phi); 
    x_value = grid_size*cos(phi);
    hypotenuse = sqrt(y_value^2 + x_value^2);

    %once we have the length we can calculate the exact coordinates it takes in
    %a 5x5 grid matrix 
    y_max_coordinate = uint8(hypotenuse*sin(theta));
    x_max_coordinate = uint8(hypotenuse*cos(theta));


    model_grid = zeros(grid_size); 

    %find which coordinate extends the furthest
    largest_vector = max(y_max_coordinate, x_max_coordinate);
    y_coordinate = y_max_coordinate;
    x_coordinate = 1; 
    for i = 1 : largest_vector
        centered_y_coordinate = y_coordinate + (grid_size - y_max_coordinate);
        model_grid(centered_y_coordinate, x_coordinate) = 1; 
        if(y_coordinate > 0) 
            y_coordinate = y_coordinate - 1; 
        end 
        if(x_coordinate < x_max_coordinate+1 ) 
            x_coordinate = x_coordinate + 1; 
        end
    end

    %once we have populated the grid with the line, we will make alterations to
    %the grid to interpret the adc value (such as with colour) 

    %initialise a RGB matrix 
    model_RGB = zeros(grid_size, grid_size, 3);

    %choose starting colour
    starting_red = 255;
    starting_green = 200; 
    starting_blue = 66;


    for colour_channel = 1 : 3
        for grid_row = 1 : grid_size
            for grid_col = 1 : grid_size
                if(model_grid(grid_row, grid_col) ~= 0)
                    switch (colour_channel)
                        case 1
                            %if we are editing the red colour channel
                            model_RGB(grid_row, grid_col, colour_channel) = starting_red;
                        case 2
                            %if we are editing the green channel, we set the
                            %green colour in proportion to the adc 
                            amount_of_green = starting_green - uint8(starting_green * normalised_adc);
                            model_RGB(grid_row, grid_col, colour_channel) = amount_of_green;
                        case 3 
                            %if we are editing the blue colour channel 
                            model_RGB(grid_row, grid_col, colour_channel) = starting_blue;
                    end
                end
            end
        end
    end
    RGB_stick_model = model_RGB;
end 