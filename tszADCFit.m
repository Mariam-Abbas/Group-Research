%Global to pass into fmincon
global data; 
global b_values; 
%Load that data into matlab from directory
data = niftiread('nonPregnant1708_2101.nii.gz'); 
protocol_21 = load('protocol_21.txt'); 

%extrapolate b-values from protocol_21 
b_values = protocol_21(:, 4); 

real_sig = getSignal(101, 91, 26, data); 
s_zero = getSZero(real_sig, b_values);

lower_bound = -1; 
upper_bound = 1; 

opt_ADC = fmincon(@optimiser, 0, [],[],[],[],lower_bound, upper_bound);

plot_signal_against(opt_ADC, 101, 91, 26);

%taking a look at the segmentation file 
[layers, coordinates] = mask_coodinates_extractor('seg_1.nii.gz');


%function to be minimised - caluclates the ssd of actual signal with
%predicted signal 

function ssd = optimiser(ADC_guess)
    global data;
    global b_values;
    
    %Obtain the original signal data of the slice we want 
    true_signal = getSignal(101, 91, 26, data);
    
    %Would like to add normalisation method right here, but not sure how to
    %normalise nii data values (Do I divide by highest value?) 
    
    %get the synthesized signal 
    synthesized_signal = getSynthSignal(true_signal, b_values, ADC_guess);
    
    ssd_temp = 0; 
    
    for i = 1 : length(true_signal)
        ssd_temp = ssd_temp + (synthesized_signal(i, 1) - true_signal(i, 1))^2;
    end
    
    ssd = ssd_temp;
end

%function to caluculate an approximation to s0 but taking mean of b = 0
%value 
function s_zero = getSZero(real_signal_data, bs)
    %First get the row numbers of when b = 0 from b_value data
    
    %We create an an empty array to store the positions of when b = 0; 
    zero_b_pos = [];
    for i = 1 : size(bs)
        if bs(i, 1) == 0
            %once we find a b_value that is equal to zero, we concatenate
            %it with the array
            zero_b_pos = [zero_b_pos, i];
        end
    end
    %now we loop through the data at the positions where b = 0 and sum up
    %the signal values; 
    sum_signal_val = 0; 
    for i = 1 : size(zero_b_pos) 
        sum_signal_val = real_signal_data(zero_b_pos(i)) + sum_signal_val;
    end
    %Then we returned the average value of the signal 
    s_zero = double(sum_signal_val / length(zero_b_pos));
end

%function to extrapolate actual data signal by passing coordinates as
%parameters 
function real_signal = getSignal(xVal, yVal, zVal, img_data)
    real_sig_temp = squeeze(img_data(xVal, yVal, zVal, :)); 
    real_signal = double(real_sig_temp); 
end

%function that produces the synthesized signal out of specified adc value
function synth_signal = getSynthSignal(og_signal, bs, adc_val)
    %get estimation of s0 to produce synthesized signal 
    s_zero = getSZero(og_signal, bs)
    %create a a zero matrix of the original signal to hold generated
    %signals
    synthesized_signal = zeros(size(og_signal));
    for i = 1 : size(bs)
        synthesized_signal(i, 1) = s_zero * exp( -1 * bs(i) * adc_val);
    end
    synth_signal = synthesized_signal;
end

%function to plot signal 
function plot_signal_against(adc_value, xVal, yVal, zVal)
    global data; 
    global b_values;
    %get the original signal 
    original_signal = getSignal(xVal, yVal, zVal, data);
    %get synthsized signal 
    synthesized_signal = getSynthSignal(original_signal, b_values, adc_value);
    %plot the graphs against each other
    figure();
    subplot(1,2,1);
    scatter(b_values, original_signal);title('original signal');
    subplot(1,2,2); 
    scatter(b_values, synthesized_signal);title('synthesized signal');
end

