function model_historgram_gen(parameter_matrix, voxel_num)
    %For each of the models we produce a histogram using the parameter
    %matrix 
    
    %set up a single display figure 
    
    %loop through all the models
    
    for model_num = 1 : 13 
       %get the BIC value matrix for all the voxel of that model 
        target_model_BICs = parameter_matrix(1:voxel_num, model_num, 13);
        
        %Then we round up the BIC values to get a nicer histogram 
        target_model_BICs = round(target_model_BICs, 0); 
        

        %afterwards we just plot the BIC values 
   
        subplot(4, 4, model_num); 
        
        histogram(target_model_BICs);
        title(model_num); 
    end 
end 
