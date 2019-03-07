%Globals:
global data;
data = niftiread("nonPregnant1708_2101.nii.gz");
realSignals = data;
trueSignal = squeeze(data(101, 91, 26, :));

global protocol_21;
protocol_21 = load('protocol_21.txt');
bValues = protocol_21(:,4);

global plotsignal;
global leftover_error;

global model_to_fit;
% 0 - Ball Model
% 1 - Ball-Ball Model
% 2 - Stick Model
% 3 - Stick-Stick Model
% 4 - IVIM Model
% 5 - Zeppelin Model
% 6 - Zeppelin-Zeppelin Model
% 7 - Ball- Stick Model
% 8 - Stick - Ball Model
% 9 - Ball - Zeppelin Model
% 10 - Zeppelin - Ball Model
% 11 - Ball - Zeppelin Model
% 12 - Zeppelin - Ball Model
% 13 - Stick - Zeppelin Model
% 14 - Zeppelin -Stick Model


% 7 - Tensor Model

%x0 = get_start_value(lb, ub);

%Run and Display Result:
    [u, coordinates] = mask_coodinates_extractor('seg_1.nii.gz');
    [x, y, z] = size(coordinates);
    for coordinate_value = 1 : 1;%x; 
        for model = 1: 1;%7;
            x_val = coordinates(1, coordinate_value);
            y_val = coordinates(2, coordinate_value);
            z_val = coordinates(3, coordinate_value);
            minimised = run_fmincon(model, x_val, y_val, z_val);
            ranking(model) = minimised;
        end 
    end
    
function minimised = run_fmincon(model_number, x, y ,z)
    global data;
    global true_signal;
    true_signal = squeeze(data(x, y, z, :));
    global model_to_fit
    switch model_number
        
        % 0 - Ball Model
        case 0
        model_to_fit = 0;
        %ADCGuess(1) = ADC ; ADCGuess(2) = S0;
        lb = [0; 0];    
        x0 = [0.01; 250];

        % 1 - Ball-Ball Model
        case 1
        model_to_fit = 1;
        %ADCGuess(1) = D ; %ADCGuess(2) = D* (Blood); ADCGuess(3) = S0; ADCGuess(4) = f
        lb = [0; 0.003; 0; 0];  
        ub = [0.003; exp(1000); exp(1000); 1];
        x0 = [0.0001; 0.005; 250; 0.5];

        % 2 - Stick Model
        case 2
        model_to_fit = 2;
        %ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = theta; ADCGuess(4) = Phi
        lb = [0; 0; -2*pi; -2*pi];  
        ub = [exp(1000); exp(1000); 2*pi; 2*pi];
        x0 = [0.003; 200; 0; 0];

        % 3 - Stick-Stick Model
        case 3
        model_to_fit = 3;
        %ADCGuess(1) = D ; %ADCGuess(2) = D* (Blood); ADCGuess(3) = S0; ADCGuess(4) = theta;
        %ADCGuess(5) = Phi; ADCGuess(6) = f 

        lb = [0; 0.003; 0; -2*pi; -2*pi; 0];  
        ub = [0.003; exp(1000); exp(1000); 2*pi; 2*pi; 1];
        x0 = [0.0001; 0.005; 250; 0; 0; 0.5];    

        % 4 - IVIM Model
        case 4
        model_to_fit = 4;
        %ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = f; ADCGuess(4) = D* (perfusioin coefficent)
        lb = [0; 0; 0; 0];  
        ub = [exp(1000); exp(1000); 1; exp(1000)];
        x0 = [0.0001; 250; 0.5; 0.01];
        
        % 5 - Zeppelin Model
        case 5
        model_to_fit = 5;
        %ADCGuess(1) = S0 ; ADCGuess(2) = theta; ADCGuess(3) = Phi ;
        %ADCGuess(4) = alpha; ADCGuess(5) = beta;
        lb = [0; -2*pi; -2*pi; 0 ; 0];  
        ub = [exp(1000); 2*pi; 2*pi; exp(1000); exp(1000)];
        x0 = [200; 0.01; 0.01; 0.01; 0.01];
        
        % 6 - Zeppelin-Zeppelin Model
        case 6
        model_to_fit = 6;
        %ADCGuess(1) = D ; %ADCGuess(2) = D* (Blood); ADCGuess(3) = S0; ADCGuess(4) = theta;
        %ADCGuess(5) = Phi; %ADCGuess(6) = alpha; ADCGuess(7) = beta;  ADCGuess(8) = theta*;
        %ADCGuess(9) = Phi*; %ADCGuess(10) = alpha*; ADCGuess(11) = beta*; ADCGuess(12) = f

        lb = [0; 0.003; 0; -2*pi; -2*pi; 0; 0; -2*pi; -2*pi; 0; 0; 0.5];  
        ub = [0.003; exp(1000); exp(1000); 2*pi; 2*pi; exp(1000); exp(1000); 2*pi; 2*pi; exp(1000); exp(1000); 1];
        x0 = [0.0001; 0.005; 250; 0; 0; 0; 0; 0; 0; 0; 0; 0.5];  
             
        % 7 - Ball-Stick Model 
        model_to_fit = 7;
        %ADCGuess(1) = D ; %ADCGuess(2) = D* (Blood); ADCGuess(3) = S0; ADCGuess(4) = theta;
        %ADCGuess(5) = Phi; ADCGuess(6) = f 
        %ADCGuess(1) = ADC ; ADCGuess(2) = S0;
        lb = [0; 0];  
        ub = [exp(1000); exp(1000)];
        x0 = [0.01; 250];


    end

    ADCopt = fmincon(@cominedOptimise, x0, [], [], [], [], lb, ub);
    disp(ADCopt);
    plot_ADC(ADCopt);
    
    global leftover_error;
    minimised = leftover_error;
end


%this is the function that gets minimised
function sum = cominedOptimise(ADCGuess)
    global plotsignal;
    
    global true_signal;
    
    global protocol_21;
    protocol_21 = double(protocol_21);
    
    global bValues;
    bValues = protocol_21(:,4);


% 7 - Ball- Stick Model
% 8 - Stick - Ball Model
% 9 - Ball - Zeppelin Model
% 10 - Zeppelin - Ball Model
% 11 - Ball - Zeppelin Model
% 12 - Zeppelin - Ball Model
% 13 - Stick - Zeppelin Model
% 14 - Zeppelin -Stick Model

    
    %Produce Test Signals
    global model_to_fit;
    switch model_to_fit
        
        % 0 - Ball Model
        case 0
            for i = 1:size(bValues)
               %ADCGuess(1) = ADC ; ADCGuess(2) = S0;
               testSignal(i) = ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1));
               plotsignal(i) = testSignal(i);
            end
            
        % 1 - Ball-Ball Model
        case 1
            %ADCGuess(1) = D ; %ADCGuess(2) = D* (Blood); ADCGuess(3) = S0; ADCGuess(4) = f
            
            for i = 1:size(bValues)              
               testSignal(i) =  ADCGuess(4)*(ADCGuess(3)*exp(-1*bValues(i)*(ADCGuess(1)))) + (1-ADCGuess(4))*(ADCGuess(3)*exp(-1*bValues(i)*ADCGuess(2)));
               plotsignal(i) = testSignal(i);
            end
        
        % 2 - Stick Model
        case 2            
            %ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = theta; ADCGuess(4) = Phi
            
            gx = protocol_21(:,1);
            gy = protocol_21(:,2);
            gz = protocol_21(:,3);
            
             for i = 1:size(bValues)
                G = [gx(i); gy(i); gz(i)];
                n = [sin(ADCGuess(3))*cos(ADCGuess(4)); sin(ADCGuess(3))* sin(ADCGuess(4)); cos(ADCGuess(3))];
                
               testSignal(i) =  ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1)*(dot(n,G))^2);
               plotsignal(i) = testSignal(i);
             end
             
       % 3 - Stick-Stick Model      
       case 3            
        %ADCGuess(1) = D ; %ADCGuess(2) = D* (Blood); ADCGuess(3) = S0; ADCGuess(4) = theta;
        %ADCGuess(5) = Phi; ADCGuess(6) = f 

        gx = protocol_21(:,1);
        gy = protocol_21(:,2);
        gz = protocol_21(:,3);

         for i = 1:size(bValues)
            G = [gx(i); gy(i); gz(i)];
            n = [sin(ADCGuess(4))*cos(ADCGuess(5)); sin(ADCGuess(4))* sin(ADCGuess(5)); cos(ADCGuess(4))];
            
            testSignal(i) =  ADCGuess(6)*(ADCGuess(3)*exp(-1*bValues(i)*ADCGuess(1))) + (1-ADCGuess(5))*(ADCGuess(3)*exp(-1*bValues(i)*ADCGuess(2)*(dot(n,G))^2));
            plotsignal(i) = testSignal(i);
         end
        
        % 4 - IVIM Model
        case 4
            %ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = f; ADCGuess(4) = D* (perfusioin coefficent)
            
            for i = 1:size(bValues)              
               testSignal(i) =  ADCGuess(3)*(ADCGuess(2)*exp(-1*bValues(i)*(ADCGuess(1) + ADCGuess(4)))) + (1-ADCGuess(3))*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1)));
               plotsignal(i) = testSignal(i);               
            end
            

        % 5 - Zeppelin Model    
        case 5           
            %ADCGuess(1) = S0 ; ADCGuess(2) = theta; ADCGuess(3) = Phi ;
            %ADCGuess(4) = alpha; ADCGuess(5) = beta;
            
            gx = protocol_21(:,1);
            gy = protocol_21(:,2);
            gz = protocol_21(:,3);
            
             for i = 1:size(bValues)
                G = [gx(i); gy(i); gz(i)];
                G_t = transpose(G);
                
                n = [sin(ADCGuess(2))*cos(ADCGuess(3)); sin(ADCGuess(2))* sin(ADCGuess(3)); cos(ADCGuess(2))];
                n_t = transpose(n);
                
              testSignal(i) =  ADCGuess(1)*exp(-1*bValues(i)*G_t*(ADCGuess(4)*n*n_t + ADCGuess(5)*eye(3))*(G));
              plotsignal(i) = testSignal(i);
             end
             
      
        % 6 - Zeppelin-Zeppelin Model
        case 6
        %ADCGuess(1) = D ; %ADCGuess(2) = D* (Blood); ADCGuess(3) = S0; ADCGuess(4) = theta;
        %ADCGuess(5) = Phi; %ADCGuess(6) = alpha; ADCGuess(7) = beta;  ADCGuess(8) = theta*;
        %ADCGuess(9) = Phi*; %ADCGuess(10) = alpha*; ADCGuess(11) = beta*; ADCGuess(12) = f
            

            gx = protocol_21(:,1);
            gy = protocol_21(:,2);
            gz = protocol_21(:,3);
            
             for i = 1:size(bValues)
                G = [gx(i); gy(i); gz(i)];
                G_t = transpose(G);
                
                n = [sin(ADCGuess(4))*cos(ADCGuess(5)); sin(ADCGuess(4))* sin(ADCGuess(5)); cos(ADCGuess(4))];
                n_t = transpose(n);
                
                D = (ADCGuess(6)*n*n_t + ADCGuess(7)*eye(3))
                
                n2 = [sin(ADCGuess(8))*cos(ADCGuess(9)); sin(ADCGuess(8))* sin(ADCGuess(9)); cos(ADCGuess(8))];
                n_t2 = transpose(n);
                
                D2 = (ADCGuess(10)*n2*n_t2 + ADCGuess(11)*eye(3))
                
              testSignal(i) =  ADCGuess(12)*(ADCGuess(3)*exp(-1*bValues(i)*G_t*D*(G))) + (1 -ADCGuess(12))*(ADCGuess(3)*exp(-1*bValues(i)*G_t*D2*(G)));
              plotsignal(i) = testSignal(i);
             end
             
       case 10000      %MRTrix TENSOR -THIS DOESNT GET USED     
            %ADCGuess(1) = S0 ; ADCGuess(2) = theta; ADCGuess(3) = Phi ;
            %ADCGuess(4) = alpha; ADCGuess(5) = beta;
            
            gx = protocol_21(:,1);
            gy = protocol_21(:,2);
            gz = protocol_21(:,3);
            
             for i = 1:size(bValues)
                G = [gx(i); gy(i); gz(i)];
                G_t = transpose(G);
                
                parallel_difussity =  (ADCGuess(4) + ADCGuess(5)); % d|| = alpha + beta
                d_up_tack1 = ADCGuess(5); %d? = beta
                d_up_tack2 = dot(parallel_difussity, d_up_tack1); %?????
                
                n = [sin(ADCGuess(2))*cos(ADCGuess(3)); sin(ADCGuess(2))* sin(ADCGuess(3)); cos(ADCGuess(2))];
                n_t = transpose(n);
                
                a_t_div_b_t = -1*(n(1)/n(2));
                a_t = a_t_div_b_t;
                b_t = 1;
                c_t = n(3);
                
                n_up_tack1 = [a_t; b_t; c_t];
                n_up_tack1_t = transpose(n_up_tack1);
                
                n_up_tack2 = cross(n, n_up_tack1);
                n_up_tack2_t = transpose(n_up_tack2);
                
                D = (parallel_difussity)*n*n_t + d_up_tack1*n_up_tack1*n_up_tack1_t + d_up_tack2*n_up_tack2*n_up_tack2_t;
                
                testSignal(i) =  ADCGuess(1)*exp(-1*bValues(i)*G_t*D*G);
                plotsignal(i) = testSignal(i);
             end
                     
    end
    
    %Find Mean Squared Error
    sum = 0;
    for i = 1:size(true_signal)
        sum = sum + (testSignal(i)-true_signal(i))^2;
    end
    sum = double(sum);
    
    global leftover_error;
    leftover_error = sum;
end

%Plots Estimated Values vs Actual Values
function foo = plot_ADC(ADCGuess)

    global plotsignal;
    global protocol_21;
    protocol_21 = double(protocol_21);    
    bValues = protocol_21(:,4);
    
    
    figure();
    hold on;
    
    global true_signal;
    
    global data; 
    scatter(bValues , plotsignal);
    scatter(bValues, true_signal);
    
end

%function that takes a mask and calculates the properties of each voxel
%ONLY FOR ADC
function ADC_map = ADC_masked_map(nifti_file, original_data)
    global data;
    global tempTrueSignal; 
    global model_to_fit;
   
    %get the coordinate of the mask
    segment_coordinates = mask_coodinates_extractor(nifti_file);
   
    %get the size of our original data
    [max_x, max_y, max_z, b] = size(original_data);
   
    %Find the layer the mask is by checking which coordinate doesn't change
    x_change = 0;
    y_change = 0;
    z_change = 0;
   
    temp_x = segment_coordinates(1, 1);
    temp_y = segment_coordinates(1, 2);
    temp_z = segment_coordinates(1, 3);
   
    %get the size of the coordinates matrix
    [i, j] = size(segment_coordinates);
    for coordinate_number = 2 : j
        new_x = segment_coordinates(1, coordinate_number);
        new_y = segment_coordinates(2, coordinate_number);
        new_z = segment_coordinates(3, coordinate_number);
       
        %check if there are changes in x, y and z
        if((x_change == 0) && (temp_x ~= new_x))
            x_change = 1;
        end
        if((y_change == 0) && (temp_y ~= new_y))
            y_change = 1;
        end
        if((z_change == 0) && (temp_z ~= new_z))
            z_change = 1;
        end
       
        %check if there is the more than one signal
        if(x_change == 1 && y_change == 1 && z_change == 1)
            disp("there is more than one layer in this segment file, abort");
            return;
        end
    end
   
    %set the plane of the segment  x = 0, y = 1, z = 2
    if(x_change == 0)
        segment_plane = 0;
    elseif(y_change == 0)
        segment_plane = 1;
    else
        segment_plane = 2;
    end
   
    %initilise a masked_canvas to map the adc onto 
    %Calculate the model fit for each voxel dependant on the plane 
    switch(segment_plane)
        case 0
            %if the segment is on the x plane
            %set the canvas of the segment mapping to be the max of the y
            %and z
            masked_canvas = zeros(max_z, max_y);
        case 1
            %if the segment is on the y plane
            %set the canvas of the segment mapping to be the max of x and z
            masked_canvas = zeros(max_x, max_z);
        case 2
            %if the segment is on the z plane
            %set the canvas of the segment mapping to be the max of x and y
            masked_canvas = zeros(max_x, max_y);
    end
   
    %Value to keep track of the largest number to be used for normalisation
    largest_adc_value = 0;
   
    %loop through each voxel from the coordinates
    for coordinate_number = 1 : j
        %get the coordinates of the coordinates
        x_coordinate = segment_coordinates(1, coordinate_number);
        y_coordinate = segment_coordinates(2, coordinate_number);
        z_coordinate = segment_coordinates(3, coordinate_number);
        %get the new true signal value
        tempTrueSignal = double(squeeze(data(x_coordinate, y_coordinate, z_coordinate, :)));
        model_to_fit = 0;
        %ADCGuess(1) = ADC ; ADCGuess(2) = S0;
        lb = [0; 0];   
        x0 = [0.01; 0.01];
        fitting_vals = fmincon(@cominedOptimise, x0, [], [], [], [], lb, []);
       
        %get the ADC
        adc_val = fitting_vals(1);
        %populate the adc value in the correct orientation
        switch (segment_plane)
            case 0
                masked_canvas(z_coordinate, y_coordinate) = adc_val;
            case 1
                masked_canvas(x_coordinate, z_coordinate) = adc_val;
            case 2
                masked_canvas(x_coordinate, y_coordinate) = adc_val;
        end
       
        %check if the new adc is larger
        if(adc_val > largest_adc_value)
            largest_adc_value = adc_val;
        end
    end
    %After looping all coordinates, we normalise the data between 0 to 1
    masked_canvas = masked_canvas / largest_adc_value;
   
    %Apply a 255 gradient mask over masked canvas
    masked_canvas = uint8(masked_canvas * 255);
    masked_canvas = imadjust(masked_canvas);
    %Output the masked canvas
    figure;
    imshow(masked_canvas);
    ADC_map = masked_canvas;
end


function x0 = get_start_value(lb, ub)

end
