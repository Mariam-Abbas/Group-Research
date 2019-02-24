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
%global best_ADC;


%Fit Models
lb = [0 0];    
x0 = [0.01; 0.01];
    
model_to_fit = 0;
ADCopt = fmincon(@cominedOptimise, x0, [], [], [], [], lb, []);

disp(ADCopt);
%disp(best_ADC);
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
    for i = 1:size(bValues)
       % testSignal(i) = s0Guess*exp(-1*bValues(i)*ADCGuess);
       testSignal(i) = ADCGuess(2,1)*exp(-1*bValues(i)*ADCGuess(1,1));
    end
    
    %Find Mean Squared Error (What is minimised0
    sum = 0;
    for i = 1:size(trueSignal)
        sum = sum + (testSignal(i)-trueSignal(i))^2;
    end
    sum = double(sum);
    
end

%function that replots the values
function foo = plot_ADC(best_ADC)
    global data; 
    global protocol_21;
    protocol_21 = double(protocol_21);    
    bValues = protocol_21(:,4);
    for i = 1:size(bValues)
       plotsignal(i) = best_ADC(2)*exp(-1*bValues(i)*best_ADC(1,1));
    end
    figure();
    hold on;
    scatter(bValues , plotsignal);
    scatter(bValues, squeeze(data(101, 91, 26, :)));
    
end
