global data;
data = niftiread("nonPregnant1708_2101.nii.gz");
global protocol_21;
realSignals = data;
bValues = protocol_21(:,4);
trueSignal = squeeze(data(101, 91, 26, :));
global best_ADC;

%testSignal = ADCFit(ADCGuess, s0Guess, bValues);
%optimise(testSignal, trueSignal);
lb = 0;    
    
ADCopt = fmincon(@cominedOptimise, [0.01], [], [], [], [], lb, [], []);
disp(ADCopt);
disp(best_ADC);
plot_ADC(best_ADC);

%dont worry about this
function testSignal = ADCFit(ADCGuess, s0Guess, bValues)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    for i = 1:size(bValues)
        testSignal(i) = s0Guess*exp(-1*bValues(i)*ADCGuess);
    end
    disp(testSignal);
end

%dont worry about this
function sum = optimise(testSignal, trueSignal)
    sum = 0;
    for i = 1:size(trueSignal)
        sum = sum + (testSignal(i)-trueSignal(i))^2;
    end
    disp(sum);
end

%this is the function that gets minimised
function sum = cominedOptimise(ADCGuess)
   % data = niftiread('nonPregnant1708_2101.nii.gz');
   global data;
    data = double(data);
    global protocol_21;
   % protocol_21 = load('protocol_21.txt');
    protocol_21 = double(protocol_21);
    trueSignal = squeeze(data(101, 91, 26, :));
    %s0Guess = 1/mean(trueSignal);
    trueSignal = trueSignal /trueSignal(1);
    global bValues;
    bValues = protocol_21(:,4);
    for i = 1:size(bValues)
       % testSignal(i) = s0Guess*exp(-1*bValues(i)*ADCGuess);
       testSignal(i) = exp(-1*bValues(i)*ADCGuess);
    end
    
    sum = 0;
    for i = 1:size(trueSignal)
        sum = sum + (testSignal(i)-trueSignal(i))^2;
    end
    sum = double(sum);
    global best_ADC;
    best_ADC = ADCGuess;
end

%function that replots the values
function foo = plot_ADC(best_ADC)
    global protocol_21;
    protocol_21 = double(protocol_21);    
    bValues = protocol_21(:,4);
    for i = 1:size(bValues)
       plotsignal(i) = exp(-1*bValues(i)*best_ADC);
    end
    figure(plot(bValues , plotsignal));
end
