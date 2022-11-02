function [x,y,ux,uy]= trajectory(A,B,tf,dt)
        coeff_x = double(coeffs(A(1),B(1),0,tf)); 
        coeff_y = double(coeffs(A(2),B(2),0,tf));
    
        num_samples = tf/dt ;
        x1 = zeros(length(num_samples),1); 
        y1 = x1; ux1 = y1; uy1 = ux1;
        
        for i = 1:num_samples 
            time = (i-1)*dt; 
            if (time < tf*0.1)%first polynom
                x1(i) = polyval(coeff_x(1:5),time);
                y1(i) = polyval(coeff_y(1:5),time);
                ux1(i) = polyval(polyder(coeff_x(1:5)),time);
                uy1(i) = polyval(polyder(coeff_y(1:5)),time);
            elseif (time < tf*0.9) % moving accoriding second polynom
                x1(i) = polyval(coeff_x(6:7),time);
                y1(i) = polyval(coeff_y(6:7),time);
                ux1(i) = polyval(polyder(coeff_x(6:7)),time);
                uy1(i) = polyval(polyder(coeff_y(6:7)),time);
            elseif (time < tf) % moving according 3rd polynom 
                x1(i) = polyval(coeff_x(8:12),time);
                y1(i) = polyval(coeff_y(8:12),time);
                ux1(i) = polyval(polyder(coeff_x(8:12)),time);
                uy1(i) = polyval(polyder(coeff_y(8:12)),time);
            end
        end
        x=x1; y=y1; ux=ux1; uy=uy1;