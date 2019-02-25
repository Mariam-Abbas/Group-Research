%Globals:
global data;
data = niftiread("nonPregnant1708_2101.nii.gz");
realSignals = data;
trueSignal = squeeze(data(101, 91, 26, :));

global protocol_21;
protocol_21 = load('protocol_21.txt');
bValues = protocol_21(:,4);

global model_to_fit;
% 0 - ADC Model
% 1 - Ball-Ball Model
% 2 - Stick Model
% 3 - Stick-Ball Model?
% 4 - IVIM Model


%Fit Models

model_to_fit = 0;
%ADCGuess(1) = ADC ; ADCGuess(2) = S0;
lb = [0; 0];    
x0 = [0.01; 0.01];
    
model_to_fit = 1;
%ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = f
lb = [0; 0; 0];  
ub = [exp(1000); exp(1000); 1];
x0 = [0.01; 0.01; 0.01];

model_to_fit = 2;
%ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = theta; ADCGuess(4) = Phi
lb = [0; 0; -2*pi; -2*pi];  
ub = [exp(1000); exp(1000); 2*pi; 2*pi];
x0 = [0.01; 200; 0.01; 0.01];

model_to_fit = 3;
%ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = theta;
%ADCGuess(4) = Phi; ADCGuess(5) = f 

lb = [0; 0; -2*pi; -2*pi; 0];  
ub = [exp(1000); exp(1000); 2*pi; 2*pi; 1];
x0 = [0.01; 200; 0.01; 0.01; 0.5];

model_to_fit = 4;
%ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = f; ADCGuess(4) = D* (perfusioin coefficent)
lb = [0; 0; 0; 0];  
ub = [exp(1000); exp(1000); 1; exp(1000)];
x0 = [0.01; 0.01; 0.01; 0.01];


%Run and Display Result:
ADCopt = fmincon(@cominedOptimise, x0, [], [], [], [], lb, ub);
disp(ADCopt);
plot_ADC(ADCopt);

%this is the function that gets minimised
function sum = cominedOptimise(ADCGuess)
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
            end

        case 1
            %ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = f
            
            for i = 1:size(bValues)              
               testSignal(i) =  ADCGuess(3)*(ADCGuess(2)*exp(-1*bValues(i)*(ADCGuess(1)))) + (1-ADCGuess(3))*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1)));
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
             end
             
       case 3            
        %ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = theta;
        %ADCGuess(4) = Phi; ADCGuess(5) = f 

        gx = protocol_21(:,1);
        gy = protocol_21(:,2);
        gz = protocol_21(:,3);

         for i = 1:size(bValues)
            G = [gx(i); gy(i); gz(i)];
            n = [sin(ADCGuess(3))*cos(ADCGuess(4)); sin(ADCGuess(3))* sin(ADCGuess(4)); cos(ADCGuess(3))];
            testSignal(i) =  ADCGuess(5)*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1))) + (1-ADCGuess(5))*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1)*(dot(n,G))^2));
         end
        
        case 4
            %ADCGuess(1) = ADC ; ADCGuess(2) = S0; ADCGuess(3) = f; ADCGuess(4) = D* (perfusioin coefficent)
            
            for i = 1:size(bValues)              
               testSignal(i) =  ADCGuess(3)*(ADCGuess(2)*exp(-1*bValues(i)*(ADCGuess(1) + ADCGuess(4)))) + (1-ADCGuess(3))*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1)));
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
    
    global protocol_21;
    protocol_21 = double(protocol_21);    
    bValues = protocol_21(:,4);
    
    global model_to_fit
    switch model_to_fit
        case 0
            for i = 1:size(bValues)
               plotsignal(i) = ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1));
            end

        case 1
            for i = 1:size(bValues)              
               plotsignal(i) =  ADCGuess(3)*(ADCGuess(2)*exp(-1*bValues(i)*(ADCGuess(1)))) + (1-ADCGuess(3))*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1)));
            end
            
        case 2
            gx = protocol_21(:,1);
            gy = protocol_21(:,2);
            gz = protocol_21(:,3);
            
             for i = 1:size(bValues)
                G = [gx(i); gy(i); gz(i)];
                n = [sin(ADCGuess(3))*cos(ADCGuess(4)); sin(ADCGuess(3))* sin(ADCGuess(4)); cos(ADCGuess(3))];
                
               plotsignal(i) =  ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1) * dot(n,G)^2);
             end

        case 3    
            gx = protocol_21(:,1);
            gy = protocol_21(:,2);
            gz = protocol_21(:,3);

             for i = 1:size(bValues)
                G = [gx(i); gy(i); gz(i)];
                n = [sin(ADCGuess(3))*cos(ADCGuess(4)); sin(ADCGuess(3))* sin(ADCGuess(4)); cos(ADCGuess(3))];
                plotsignal(i) =  ADCGuess(5)*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1))) + (1-ADCGuess(5))*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1)*(dot(n,G))^2));
             end
        
        case 4
           for i = 1:size(bValues)              
               plotsignal(i) =  ADCGuess(3)*(ADCGuess(2)*exp(-1*bValues(i)*(ADCGuess(1) + ADCGuess(4)))) + (1-ADCGuess(3))*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1)));
           end
        
    end
    figure();
    hold on;
    
    global data; 
    scatter(bValues , plotsignal);
    scatter(bValues, squeeze(data(101, 91, 26, :)));
    
end
