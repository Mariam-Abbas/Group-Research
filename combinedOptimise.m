%this is the function that gets minimised

function sum_1 = combinedOptimise(ADCGuess)
    global plotsignal;
    
    global protocol_21;
    protocol_21 = double(protocol_21);
    
    bValues = protocol_21(:,4);
    
    testSignal = 0;

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
        %ADCGuess(5) = Phi; ADCGuess(6) = theta; %ADCGuess(7) = Phi ADCGuess(8) = f 

        gx = protocol_21(:,1);
        gy = protocol_21(:,2);
        gz = protocol_21(:,3);

         for i = 1:size(bValues)
            G = [gx(i); gy(i); gz(i)];
            
            n1 = [sin(ADCGuess(4))*cos(ADCGuess(5)); sin(ADCGuess(4))* sin(ADCGuess(5)); cos(ADCGuess(4))];
            n2 = [sin(ADCGuess(6))*cos(ADCGuess(7)); sin(ADCGuess(6))* sin(ADCGuess(7)); cos(ADCGuess(6))];
            
            testSignal(i) =  ADCGuess(8)*(ADCGuess(3)*exp(-1*bValues(i)*ADCGuess(1)*(dot(n1,G))^2)) + (1-ADCGuess(8))*(ADCGuess(3)*exp(-1*bValues(i)*ADCGuess(2)*(dot(n2,G))^2));
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
             
        % 7 - Ball-Stick Model 
        case 7
        %ADCGuess(1) = D ; %ADCGuess(2) = D* (Blood); ADCGuess(3) = S0; ADCGuess(4) = theta;
        %ADCGuess(5) = Phi; ADCGuess(6) = f ;
            gx = protocol_21(:,1);
            gy = protocol_21(:,2);
            gz = protocol_21(:,3);

             for i = 1:size(bValues)
                G = [gx(i); gy(i); gz(i)];
                n = [sin(ADCGuess(4))*cos(ADCGuess(5)); sin(ADCGuess(4))* sin(ADCGuess(5)); cos(ADCGuess(4))];

                testSignal(i) =  ADCGuess(6)*(ADCGuess(3)*exp(-1*bValues(i)*ADCGuess(1))) + (1-ADCGuess(6))*(ADCGuess(3)*exp(-1*bValues(i)*ADCGuess(2)*(dot(n,G))^2));
                plotsignal(i) = testSignal(i);
             end
        
  
        % 8 - Stick-Ball Model   
        case 8
            %ADCGuess(1) = D ; %ADCGuess(2) = D* (Blood); ADCGuess(3) = S0; ADCGuess(4) = theta;
            %ADCGuess(5) = Phi; ADCGuess(6) = f ;
            gx = protocol_21(:,1);
            gy = protocol_21(:,2);
            gz = protocol_21(:,3);

             for i = 1:size(bValues)
                G = [gx(i); gy(i); gz(i)];
                n = [sin(ADCGuess(4))*cos(ADCGuess(5)); sin(ADCGuess(4))* sin(ADCGuess(5)); cos(ADCGuess(4))];

                testSignal(i) =  ADCGuess(6)*(ADCGuess(3)*exp(-1*bValues(i)*ADCGuess(1)*(dot(n,G))^2)) + (1-ADCGuess(6))*(ADCGuess(3)*exp(-1*bValues(i)*ADCGuess(2)));
                plotsignal(i) = testSignal(i);
             end
             
        % 9 - Ball-Zeppelin Model
        case 9
        %ADCGuess(1) = D ; ADCGuess(2) = S0; ADCGuess(3) = theta;
        %ADCGuess(4) = Phi; %ADCGuess(5) = alpha; ADCGuess(6) = beta; ADCGuess(7) = f
        
            gx = protocol_21(:,1);
            gy = protocol_21(:,2);
            gz = protocol_21(:,3);
            
             for i = 1:size(bValues)
                G = [gx(i); gy(i); gz(i)];
                G_t = transpose(G);
                
                n = [sin(ADCGuess(3))*cos(ADCGuess(4)); sin(ADCGuess(3))* sin(ADCGuess(4)); cos(ADCGuess(3))];
                n_t = transpose(n);
                
                D = (ADCGuess(5)*n*n_t + ADCGuess(6)*eye(3));
                
              testSignal(i) =  ADCGuess(7)*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1))) + (1 -ADCGuess(7))*(ADCGuess(2)*exp(-1*bValues(i)*G_t*D*(G)));
              plotsignal(i) = testSignal(i);
             end
             
        % 10 - Zeppelin-Ball Model
        case 10
        %ADCGuess(1) = D* (Blood); ADCGuess(2) = S0; ADCGuess(3) = theta;
        %ADCGuess(4) = Phi; %ADCGuess(5) = alpha; ADCGuess(6) = beta; ADCGuess(7) = f
        
            gx = protocol_21(:,1);
            gy = protocol_21(:,2);
            gz = protocol_21(:,3);
            
             for i = 1:size(bValues)
                G = [gx(i); gy(i); gz(i)];
                G_t = transpose(G);
                
                n = [sin(ADCGuess(3))*cos(ADCGuess(4)); sin(ADCGuess(3))* sin(ADCGuess(4)); cos(ADCGuess(3))];
                n_t = transpose(n);
                
                D = (ADCGuess(5)*n*n_t + ADCGuess(6)*eye(3));
                
              testSignal(i) =  ADCGuess(7)*(ADCGuess(2)*exp(-1*bValues(i)*G_t*D*(G))) + (1 -ADCGuess(7))*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1)));
              plotsignal(i) = testSignal(i);
             end

      % 11 - Stick-Zeppelin Model 
       case 11  
       %ADCGuess(1) = D ; ADCGuess(2) = S0; ADCGuess(3) = theta; %ADCGuess(4) = Phi; %ADCGuess(5) = alpha; 
       %ADCGuess(6) = beta; ADCGuess(7) = theta; ADCGuess(8) = Phi; ADCGuess(9) = f    
           gx = protocol_21(:,1);
            gy = protocol_21(:,2);
            gz = protocol_21(:,3);
            
             for i = 1:size(bValues)
                G = [gx(i); gy(i); gz(i)];
                G_t = transpose(G);
                
                n = [sin(ADCGuess(3))*cos(ADCGuess(4)); sin(ADCGuess(3))* sin(ADCGuess(4)); cos(ADCGuess(3))];
                n_t = transpose(n);
                
                n2 = [sin(ADCGuess(7))*cos(ADCGuess(8)); sin(ADCGuess(7))* sin(ADCGuess(8)); cos(ADCGuess(7))];
                
                D = (ADCGuess(5)*n*n_t + ADCGuess(6)*eye(3));
                
              testSignal(i) =  ADCGuess(9)*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1)*(dot(n2,G))^2)) + (1 -ADCGuess(9))*(ADCGuess(2)*exp(-1*bValues(i)*G_t*D*(G)));
              plotsignal(i) = testSignal(i);
             end
             
             
        % 11 - Zeppelin-Stick Model      
        case 12    
        %ADCGuess(1) = D* (Blood) ; ADCGuess(2) = S0; ADCGuess(3) = theta;%ADCGuess(4) = Phi; %ADCGuess(5) = alpha; 
        %ADCGuess(6) = beta; ADCGuess(7) = theta; ADCGuess(8) = Phi; ADCGuess(9) = f      
            gx = protocol_21(:,1);
            gy = protocol_21(:,2);
            gz = protocol_21(:,3);
            
             for i = 1:size(bValues)
                G = [gx(i); gy(i); gz(i)];
                G_t = transpose(G);
                
                n = [sin(ADCGuess(3))*cos(ADCGuess(4)); sin(ADCGuess(3))* sin(ADCGuess(4)); cos(ADCGuess(3))];
                n_t = transpose(n);
                
                n2 = [sin(ADCGuess(7))*cos(ADCGuess(8)); sin(ADCGuess(7))* sin(ADCGuess(8)); cos(ADCGuess(7))];
                
                D = (ADCGuess(5)*n*n_t + ADCGuess(6)*eye(3));
                
              testSignal(i) =  ADCGuess(9)*(ADCGuess(2)*exp(-1*bValues(i)*G_t*D*(G))) + (1 -ADCGuess(9))*(ADCGuess(2)*exp(-1*bValues(i)*ADCGuess(1)*(dot(n2,G))^2));
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
    global true_signal;
    sum_1 = double(0);
    for i = 1:size(true_signal)
        dif_squared = double((testSignal(i)-true_signal(i))^2);
        sum_1 = sum_1 + dif_squared;
    end
    
    sum_1 = double(sum_1/51);  
    
    global leftover_error;
    leftover_error = sum_1;
    
    sum_1 = double(sum_1);  
    
    figure; 
    plot(bValues, true_signal, 'o'); 
    hold on; 
    plot(bValues, testSignal, 'o');

end
