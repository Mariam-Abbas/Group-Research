%Globals:
global data;
data = niftiread("nonPregnant1708_2101.nii.gz");
data = double(data);

global protocol_21;
protocol_21 = load('protocol_21.txt');


global plotsignal;

global leftover_error;
leftover_error = int64(0);

global model_to_fit;

global BIC;
global param; 


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
% 11 - Stick - Zeppelin Model
% 12- Zeppelin -Stick Model


% 1000 - Tensor Model %NOT DONE 
%Run and Display Result:
    [u, coordinates] = mask_coodinates_extractor('seg_1.nii.gz');  
    
    %3D array that contains 1d: Voxel, 2d: Model Type, 3d: Voxel
    %properties(13, to store max degrees of freedom and minimised + 3more for coordinates)
    
    sigma = 1/20; 
    ranking = zeros(1, 13);
    number_of_voxels = 1;
       for number_of_compartments = 3: 5
           %parameter_map{model+1} = zeros(max_y, max_x, max_z, 20); 
           for coordinate_value = 1 : number_of_voxels;
     
            x_val = coordinates(1, coordinate_value);
            y_val = coordinates(2, coordinate_value);
            z_val = coordinates(3, coordinate_value);
            minimised = int64(0);
            %[minimised, returned_parameters] = run_fmincon(model, 101, 101, 36);% NICE VOXHAL FOR NOW x_val, y_val, z_val);
            disp(coordinate_value);
            %returned_parameters = zeros(12, 1);
            [minimised, returned_parameters] = run_fmincon(number_of_compartments, x_val, y_val, z_val);
            
            BIC = log(51)*length(returned_parameters) - 2*(-1/(2*(sigma)^2))*minimised;

           end
       end

    
    %plot the historgram models 
    %model_historgram_gen(parameter_matrix, number_of_voxels);
    
    
function [minimised, fitted_parameters] = run_fmincon(number_of_compartments, x, y ,z)
    global data;
    
    global true_signal;
    true_signal = squeeze(data(x, y, z, :));
    
    global compartments;
    compartments = number_of_compartments;
    global param; 
    
    switch number_of_compartments
        case 1
            %ADCGuess(1) = S0 ; ADCGuess(2) = theta; ADCGuess(3) = Phi ;
            %ADCGuess(4) = alpha; ADCGuess(5) = beta;
            lb = [0; -2*pi; -2*pi; 0 ; 0];  
            ub = [exp(1000); 2*pi; 2*pi; exp(1000); exp(1000)];
            x0 = [true_signal(1); 1; -1; 0.003; 0.003]; % Zeppelin 5
            
        case 2
        %ADCGuess(1) = S0; ADCGuess(2) = theta; %ADCGuess(3) = Phi; %ADCGuess(4) = alpha; ADCGuess(5) = beta;  
        %ADCGuess(6) = theta*; %ADCGuess(7) = Phi*; %ADCGuess(8) = alpha*; ADCGuess(9) = beta*; ADCGuess(10) = f

        lb = [0; -2*pi; -2*pi; 0; 0; -2*pi; -2*pi; 0; 0; 0.5];  
        ub = [exp(1000); 2*pi; 2*pi; exp(1000); exp(1000); 2*pi; 2*pi; exp(1000); exp(1000); 1];
        x0 = [true_signal(1); 1; -1; 0.003; 0.003; 1; -1; 0.003; 0.003; 0.5];  %zeppelin-zeppein 6  
        
        case 3
            %ADCGuess(1) = S0; ADCGuess(2) = theta; %ADCGuess(3) = Phi; %ADCGuess(4) = alpha; ADCGuess(5) = beta;  
            %ADCGuess(6) = theta*; %ADCGuess(7) = Phi*; %ADCGuess(8) = %alpha*; ADCGuess(9) = beta*; 
            %ADCGuess(10) = theta3; %ADCGuess(11) = Phi3; %ADCGuess(12) = %alpha3; ADCGuess(13) = beta3;
            %ADCGuess(14) = f1; %ADCGuess(15) = f2; ADCGuess(16) = f3;

            lb = [0; -2*pi; -2*pi; 0; 0; -2*pi; -2*pi; 0; 0;-2*pi; -2*pi; 0; 0; 0; 0; 0];  
            ub = [exp(1000); 2*pi; 2*pi; exp(1000); exp(1000); 2*pi; 2*pi; exp(1000); exp(1000);  2*pi; 2*pi; exp(1000); exp(1000); 100;100;100];
            x0 = [true_signal(1); 1; -1; 0.003; 0.003; 1; -1; 0.003; 0.003; 1; -1; 0.003; 0.003; 50; 50; 50];  
        
        case 4           
            %ADCGuess(1) = S0; ADCGuess(2) = theta; %ADCGuess(3) = Phi; %ADCGuess(4) = alpha; ADCGuess(5) = beta;  
            %ADCGuess(6) = theta*; %ADCGuess(7) = Phi*; %ADCGuess(8) = %alpha*; ADCGuess(9) = beta*; 
            %ADCGuess(10) = theta3; %ADCGuess(11) = Phi3; %ADCGuess(12) = %alpha3; ADCGuess(13) = beta3;
            %ADCGuess(14) = theta4; %ADCGuess(15) = Phi4; %ADCGuess(16) = %alpha4; ADCGuess(17) = beta4;
            %ADCGuess(18) = f1; %ADCGuess(19) = f2; ADCGuess(20) = f3;
            %ADCGuess(21) = f4

            lb = [0; -2*pi; -2*pi; 0; 0; -2*pi; -2*pi; 0; 0;-2*pi; -2*pi; 0; 0; -2*pi; -2*pi; 0; 0; 0; 0; 0; 0];  
            ub = [exp(1000); 2*pi; 2*pi; exp(1000); exp(1000); 2*pi; 2*pi; exp(1000); exp(1000);  2*pi; 2*pi; exp(1000); exp(1000);2*pi; 2*pi; exp(1000); exp(1000); 100;100;100; 100];
            x0 = [true_signal(1); 1; -1; 0.003; 0.003; 1; -1; 0.003; 0.003; 1; -1; 0.003; 0.003; 1; -1; 0.003; 0.003; 50; 50; 50; 50]; 
        
        case 5           
            %ADCGuess(1) = S0; ADCGuess(2) = theta; %ADCGuess(3) = Phi; %ADCGuess(4) = alpha; ADCGuess(5) = beta;  
            %ADCGuess(6) = theta*; %ADCGuess(7) = Phi*; %ADCGuess(8) = %alpha*; ADCGuess(9) = beta*; 
            %ADCGuess(10) = theta3; %ADCGuess(11) = Phi3; %ADCGuess(12) = %alpha3; ADCGuess(13) = beta3;
            %ADCGuess(14) = theta4; %ADCGuess(15) = Phi4; %ADCGuess(16) = %alpha4; ADCGuess(17) = beta4;
            %ADCGuess(18) = theta5; %ADCGuess(19) = Phi5; %ADCGuess(20) = %alpha5; ADCGuess(21) = beta5;
            %ADCGuess(22) = f1; %ADCGuess(23) = f2; ADCGuess(24) = f3;
            %ADCGuess(25) = f4; %ADCGuess(26) = f5

            lb = [0; -2*pi; -2*pi; 0; 0; -2*pi; -2*pi; 0; 0;-2*pi; -2*pi; 0; 0; -2*pi; -2*pi; 0; 0; 2*pi; -2*pi; 0; 0; 0; 0; 0; 0; 0];  
            ub = [exp(1000); 2*pi; 2*pi; exp(1000); exp(1000); 2*pi; 2*pi; exp(1000); exp(1000); 2*pi; 2*pi; exp(1000); exp(1000);  2*pi; 2*pi; exp(1000); exp(1000);2*pi; 2*pi; exp(1000); exp(1000); 100;100;100;100;100];
            x0 = [true_signal(1); 1; -1; 0.003; 0.003; 1; -1; 0.003; 0.003; 1; -1; 0.003; 0.003; 1; -1; 0.003; 0.003; 1; -1; 0.003; 0.003; 50; 50; 50; 50; 50];             
    end
    options = optimset('Display','off');
    %options = struct('MaxFunEvals', 5000);
    for i = 1: 3 
        ADCopt = fmincon(@cominedOptimise, x0, [], [], [], [], lb, ub, [], options);
        x0 = param;
    end 
    fitted_parameters = param;
    %disp(ADCopt);
    %plot_ADC();
    
    global leftover_error;
    minimised = leftover_error;
end


%this is the function that gets minimised
function sum_1 = cominedOptimise(ADCGuess)
    global plotsignal;
    global protocol_21;
    protocol_21 = double(protocol_21);
    
    bValues = protocol_21(:,4);
    
    testSignal = 0;
    global compartments;
    switch compartments
        case 1
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
            
        case 2
        %ADCGuess(1) = S0; ADCGuess(2) = theta; %ADCGuess(3) = Phi; %ADCGuess(4) = alpha; ADCGuess(5) = beta;  
        %ADCGuess(6) = theta*; %ADCGuess(7) = Phi*; %ADCGuess(8) = alpha*; ADCGuess(9) = beta*; ADCGuess(10) = f

       
    gx = protocol_21(:,1);
    gy = protocol_21(:,2);
    gz = protocol_21(:,3);

     for i = 1:size(bValues)
        G = [gx(i); gy(i); gz(i)];
        G_t = transpose(G);

        n = [sin(ADCGuess(2))*cos(ADCGuess(3)); sin(ADCGuess(2))* sin(ADCGuess(3)); cos(ADCGuess(2))];
        n_t = transpose(n);

        D = (ADCGuess(4)*n*n_t + ADCGuess(5)*eye(3));

        n2 = [sin(ADCGuess(6))*cos(ADCGuess(7)); sin(ADCGuess(6))* sin(ADCGuess(7)); cos(ADCGuess(6))];
        n_t2 = transpose(n);

        D2 = (ADCGuess(8)*n2*n_t2 + ADCGuess(9)*eye(3));

      testSignal(i) =  ADCGuess(10)*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D*(G))) + (1 -ADCGuess(10))*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D2*(G)));
      plotsignal(i) = testSignal(i);
     end
        
        case 3
            %ADCGuess(1) = S0; ADCGuess(2) = theta; %ADCGuess(3) = Phi; %ADCGuess(4) = alpha; ADCGuess(5) = beta;  
            
  

    gx = protocol_21(:,1);
    gy = protocol_21(:,2);
    gz = protocol_21(:,3);

     for i = 1:size(bValues)
        G = [gx(i); gy(i); gz(i)];
        G_t = transpose(G);

        n = [sin(ADCGuess(2))*cos(ADCGuess(3)); sin(ADCGuess(2))* sin(ADCGuess(3)); cos(ADCGuess(2))];
        n_t = transpose(n);

        D = (ADCGuess(4)*n*n_t + ADCGuess(5)*eye(3));

        n2 = [sin(ADCGuess(6))*cos(ADCGuess(7)); sin(ADCGuess(6))* sin(ADCGuess(7)); cos(ADCGuess(6))];
        n_t2 = transpose(n);

        D2 = (ADCGuess(8)*n2*n_t2 + ADCGuess(9)*eye(3));
        
        %ADCGuess(6) = theta*; %ADCGuess(7) = Phi*; %ADCGuess(8) = %alpha*; ADCGuess(9) = beta*; 
          %ADCGuess(10) = theta3; %ADCGuess(11) = Phi3; %ADCGuess(12) = %alpha3; ADCGuess(13) = beta3;
         %ADCGuess(14) = f1; %ADCGuess(15) = f2; ADCGuess(16) = f3;
        
        n3 = [sin(ADCGuess(10))*cos(ADCGuess(11)); sin(ADCGuess(10))* sin(ADCGuess(11)); cos(ADCGuess(10))];
        n_t3 = transpose(n);

        D3 = (ADCGuess(12)*n3*n_t3 + ADCGuess(13)*eye(3));
        
        sum = ADCGuess(14) + ADCGuess(15) + ADCGuess(16);
        f1 = ADCGuess(14)/sum;
        f2 = ADCGuess(15)/sum;
        f3 = ADCGuess(16)/sum;
        
        testSignal(i) =  f1*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D*(G))) + f2*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D2*(G))) + f3*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D3*(G)));
        plotsignal(i) = testSignal(i);
     end
        
        case 4           
                %ADCGuess(1) = S0; ADCGuess(2) = theta; %ADCGuess(3) = Phi; %ADCGuess(4) = alpha; ADCGuess(5) = beta;  

        gx = protocol_21(:,1);
        gy = protocol_21(:,2);
        gz = protocol_21(:,3);

         for i = 1:size(bValues)
            G = [gx(i); gy(i); gz(i)];
            G_t = transpose(G);

            n = [sin(ADCGuess(2))*cos(ADCGuess(3)); sin(ADCGuess(2))* sin(ADCGuess(3)); cos(ADCGuess(2))];
            n_t = transpose(n);

            D = (ADCGuess(4)*n*n_t + ADCGuess(5)*eye(3));

            n2 = [sin(ADCGuess(6))*cos(ADCGuess(7)); sin(ADCGuess(6))* sin(ADCGuess(7)); cos(ADCGuess(6))];
            n_t2 = transpose(n);

            D2 = (ADCGuess(8)*n2*n_t2 + ADCGuess(9)*eye(3));

                %ADCGuess(6) = theta*; %ADCGuess(7) = Phi*; %ADCGuess(8) = %alpha*; ADCGuess(9) = beta*; 
            %ADCGuess(10) = theta3; %ADCGuess(11) = Phi3; %ADCGuess(12) = %alpha3; ADCGuess(13) = beta3;
            %ADCGuess(14) = theta4; %ADCGuess(15) = Phi4; %ADCGuess(16) = %alpha4; ADCGuess(17) = beta4;
            %ADCGuess(18) = f1; %ADCGuess(19) = f2; ADCGuess(20) = f3;
            %ADCGuess(21) = f4

            n3 = [sin(ADCGuess(10))*cos(ADCGuess(11)); sin(ADCGuess(10))* sin(ADCGuess(11)); cos(ADCGuess(10))];
            n_t3 = transpose(n);

            D3 = (ADCGuess(12)*n3*n_t3 + ADCGuess(13)*eye(3));
 
            n4 = [sin(ADCGuess(14))*cos(ADCGuess(15)); sin(ADCGuess(14))* sin(ADCGuess(15)); cos(ADCGuess(14))];
            n_t4 = transpose(n);
            D4 = (ADCGuess(16)*n4*n_t4 + ADCGuess(17)*eye(3));

            sum = ADCGuess(18) + ADCGuess(19) + ADCGuess(20) + ADCGuess(21);
            f1 = ADCGuess(18)/sum;
            f2 = ADCGuess(19)/sum;
            f3 = ADCGuess(20)/sum;
            f4 = ADCGuess(21)/sum;

            testSignal(i) =  f1*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D*(G))) + f2*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D2*(G))) + f3*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D3*(G))) + f4*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D4*(G)));
              plotsignal(i) = testSignal(i);
         end
        case 4           
                %ADCGuess(1) = S0; ADCGuess(2) = theta; %ADCGuess(3) = Phi; %ADCGuess(4) = alpha; ADCGuess(5) = beta;  

        gx = protocol_21(:,1);
        gy = protocol_21(:,2);
        gz = protocol_21(:,3);

         for i = 1:size(bValues)
            G = [gx(i); gy(i); gz(i)];
            G_t = transpose(G);

            n = [sin(ADCGuess(2))*cos(ADCGuess(3)); sin(ADCGuess(2))* sin(ADCGuess(3)); cos(ADCGuess(2))];
            n_t = transpose(n);

            D = (ADCGuess(4)*n*n_t + ADCGuess(5)*eye(3));

            n2 = [sin(ADCGuess(6))*cos(ADCGuess(7)); sin(ADCGuess(6))* sin(ADCGuess(7)); cos(ADCGuess(6))];
            n_t2 = transpose(n);

            D2 = (ADCGuess(8)*n2*n_t2 + ADCGuess(9)*eye(3));

            %ADCGuess(10) = theta3; %ADCGuess(11) = Phi3; %ADCGuess(12) = %alpha3; ADCGuess(13) = beta3;
            %ADCGuess(14) = theta4; %ADCGuess(15) = Phi4; %ADCGuess(16) = %alpha4; ADCGuess(17) = beta4;
            %ADCGuess(18) = theta5; %ADCGuess(19) = Phi5; %ADCGuess(20) = %alpha5; ADCGuess(21) = beta5;
            %ADCGuess(22) = f1; %ADCGuess(23) = f2; ADCGuess(24) = f3;
            %ADCGuess(25) = f4; %ADCGuess(26) = f5

            n3 = [sin(ADCGuess(10))*cos(ADCGuess(11)); sin(ADCGuess(10))* sin(ADCGuess(11)); cos(ADCGuess(10))];
            n_t3 = transpose(n);

            D3 = (ADCGuess(12)*n3*n_t3 + ADCGuess(13)*eye(3));
 
            n4 = [sin(ADCGuess(14))*cos(ADCGuess(15)); sin(ADCGuess(14))* sin(ADCGuess(15)); cos(ADCGuess(14))];
            n_t4 = transpose(n);
            D4 = (ADCGuess(16)*n4*n_t4 + ADCGuess(17)*eye(3));
            
            n5 = [sin(ADCGuess(18))*cos(ADCGuess(19)); sin(ADCGuess(18))* sin(ADCGuess(19)); cos(ADCGuess(18))];
            n_t5 = transpose(n);
            D5 = (ADCGuess(20)*n5*n_t5 + ADCGuess(21)*eye(3));

            sum = ADCGuess(22) + ADCGuess(23) + ADCGuess(24) + ADCGuess(25) + ADCGuess(26);
            f1 = ADCGuess(22)/sum;
            f2 = ADCGuess(23)/sum;
            f3 = ADCGuess(24)/sum;
            f4 = ADCGuess(25)/sum;
            f5 = ADCGuess(26)/sum;
            
            testSignal(i) =  f1*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D*(G))) + f2*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D2*(G))) + f3*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D3*(G))) + f4*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D4*(G)))+ f5*(ADCGuess(1)*exp(-1*bValues(i)*G_t*D5*(G)));   
    
         end
    end
%ADCGuess(1) = S0; ADCGuess(2) = theta; %ADCGuess(3) = Phi; %ADCGuess(4) = alpha; ADCGuess(5) = beta;  
%ADCGuess(6) = theta*; %ADCGuess(7) = Phi*; %ADCGuess(8) = alpha*; ADCGuess(9) = beta*; ADCGuess(10) = f            

     
      
    
    %Find Mean Squared Error
    global true_signal;
    sum_1 = double(0);
    for i = 1:size(true_signal)
        dif_squared = double((testSignal(i)-true_signal(i))^2);
        sum_1 = sum_1 + dif_squared;
    end
    
    global param; 
    param = ADCGuess;
    global leftover_error;
    leftover_error = sum_1;
    
    sum_1 = double(sum_1);  

end

%Plots Estimated Values vs Actual Values
function foo = plot_ADC()

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




function x0 = get_start_value(lb, ub)

end