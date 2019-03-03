%Globals:
global data;
data = niftiread("nonPregnant1708_2101.nii.gz");
realSignals = data;
trueSignal = squeeze(data(101, 91, 26, :));

global protocol_21;
protocol_21 = load('protocol_21.txt');
bValues = protocol_21(:,4);

global plotsignal;

global model_to_fit;
% 0 - ADC Model
% 1 - ADC-ADC Model?
% 2 - Stick Model
% 3 - Stick-Stick Model
% 4 - IVIM Model
% 5 - Ball Model
% 6 - Zeppelin Model
% 7 - Tensor Model

%Fit Models

model_to_fit = 0;
%ADCGuess(1) = ADC ; ADCGuess(2) = S0;
lb = [0; 0];    
x0 = [0.01; 250];

model_to_fit = 1;
%ADCGuess(1) = D ; %ADCGuess(2) = D* (Blood); ADCGuess(3) = S0; ADCGuess(4) = f
lb = [0; 0.003; 0; 0];  
ub = [0.003; exp(1000); exp(1000); 1];
x0 = [0.0001; 0.005; 250; 0.5];

model_to_fit = 2;
%ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = theta; ADCGuess(4) = Phi
lb = [0; 0; -2*pi; -2*pi];  
ub = [exp(1000); exp(1000); 2*pi; 2*pi];
x0 = [0.003; 200; 0; 0];

model_to_fit = 3;
%ADCGuess(1) = D ; %ADCGuess(2) = D* (Blood); ADCGuess(3) = S0; ADCGuess(4) = theta;
%ADCGuess(5) = Phi; ADCGuess(6) = f 

lb = [0; 0.003; 0; -2*pi; -2*pi; 0];  
ub = [0.003; exp(1000); exp(1000); 2*pi; 2*pi; 1];
x0 = [0.0001; 0.005; 250; 0; 0; 0.5];    

model_to_fit = 4;
%ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = f; ADCGuess(4) = D* (perfusioin coefficent)
lb = [0; 0; 0; 0];  
ub = [exp(1000); exp(1000); 1; exp(1000)];
x0 = [0.0001; 250; 0.5; 0.01];

model_to_fit = 5;
%ADCGuess(1) = ADC ; ADCGuess(2) = S0;
lb = [0; 0];  
ub = [exp(1000); exp(1000)];
x0 = [0.01; 250];

model_to_fit = 6;
%ADCGuess(1) = S0 ; ADCGuess(2) = theta; ADCGuess(3) = Phi ;
%ADCGuess(4) = alpha; ADCGuess(5) = beta;
lb = [0; -2*pi; -2*pi; 0 ; 0];  
ub = [exp(1000); 2*pi; 2*pi; exp(1000); exp(1000)];
x0 = [200; 0.01; 0.01; 0.01; 0.01];

model_to_fit = 7;
%ADCGuess(1) = S0 ; ADCGuess(2) = theta; ADCGuess(3) = Phi ;
%ADCGuess(4) = alpha; ADCGuess(5) = beta;
lb = [0; -2*pi; -2*pi; 0 ; 0];  
ub = [exp(1000); 2*pi; 2*pi; exp(1000); exp(1000)];
x0 = [200; 0.01; 0.01; 0.01; 0.01];



%x0 = get_start_value(lb, ub);

%Run and Display Result:
ADCopt = fmincon(@cominedOptimise, x0, [], [], [], [], lb, ub);
disp(ADCopt);
plot_ADC(ADCopt);

%this is the function that gets minimised
function sum = cominedOptimise(ADCGuess)
    global plotsignal;
    
    global data;
    data = double(data);
    
    global protocol_21;
    protocol_21 = double(protocol_21);
    trueSignal = squeeze(data(101, 91, 26, :));
    
    global bValues;
    bValues = protocol_21(:,4);
    
    %Produce Test Signals
    global model_to_fit;
    switch model_to_fit
        case 0
            for i = 1:size(bValues)
               %ADCGuess(1) = ADC ; ADCGuess(2) = S0;
               testSignal(i) = ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1));
               plotsignal(i) = testSignal(i);
            end

        case 1
            %ADCGuess(1) = D ; %ADCGuess(2) = D* (Blood); ADCGuess(3) = S0; ADCGuess(4) = f
            
            for i = 1:size(bValues)              
               testSignal(i) =  ADCGuess(4)*(ADCGuess(3)*exp(-1*bValues(i)*(ADCGuess(1)))) + (1-ADCGuess(4))*(ADCGuess(3)*exp(-1*bValues(i)*ADCGuess(2)));
               plotsignal(i) = testSignal(i);
            end
            
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
        
        case 4
            %ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = f; ADCGuess(4) = D* (perfusioin coefficent)
            
            for i = 1:size(bValues)              
               testSignal(i) =  ADCGuess(3)*(ADCGuess(2)*exp(-1*bValues(i)*(ADCGuess(1) + ADCGuess(4)))) + (1-ADCGuess(3))*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1)));
               plotsignal(i) = testSignal(i);               
            end
            
          case 5            
            %ADCGuess(1) = ADC ; ADCGuess(2) = S0; 
            
            gx = protocol_21(:,1);
            gy = protocol_21(:,2);
            gz = protocol_21(:,3);
            
             for i = 1:size(bValues)
                G = [gx(i); gy(i); gz(i)];
                G_t = transpose(G);
                
               testSignal(i) =  ADCGuess(2)*exp(-1*bValues(i)*G_t*ADCGuess(1)*eye(3)*(G));
               plotsignal(i) = testSignal(i);
             end
             
        case 6           
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
             
        case 7           
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
    for i = 1:size(trueSignal)
        sum = sum + (testSignal(i)-trueSignal(i))^2;
    end
    sum = double(sum);
    
end

%Plots Estimated Values vs Actual Values
function foo = plot_ADC(ADCGuess)

    global plotsignal;
    global protocol_21;
    protocol_21 = double(protocol_21);    
    bValues = protocol_21(:,4);
    
    
    figure();
    hold on;
    
    global data; 
    scatter(bValues , plotsignal);
    scatter(bValues, squeeze(data(101, 91, 26, :)));
    
end

function x0 = get_start_value(lb, ub)

end
