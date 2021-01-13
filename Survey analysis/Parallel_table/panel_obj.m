%% Objective function for each combination of parameter values under consideration
function F = panel_obj(x,respondent_num,table,R_type,weight_type,paramhat,cdf,reg,lb,ub)
    
    walk_time = -paramhat(1);
    wait_time = -paramhat(2);
    transit_time = -paramhat(3);
    exc_time = -paramhat(4);
    pool_time = -paramhat(5);
    price = -paramhat(6);
    ASC_transit = 0; % Baseline
    ASC_exclusive = paramhat(7);
    ASC_pooled = paramhat(8);
    
    coeff_vec = [walk_time wait_time transit_time exc_time pool_time price ASC_transit ASC_exclusive ASC_pooled];
        
    % Normalize lambda so that all 4 parameters lie between 0 and 1
    x(5) = x(5) * 100;
       
    ref = @(u1,u2,p,CE,i) reference(u1,u2,p,R_type,CE,table,i,paramhat);
    calc_error = @(a,b,c,d,e,CE_sR,u1,u2,p,R) error_func(a,b,c,d,e,CE_sR,u1,u2,p,R,weight_type,cdf);
    
    i = respondent_num;
    transit_walk = table.transit_walk_1(i) + table.transit_walk_2(i);
    transit_cost = table.Transit_cost(i);

    distance = table.Distance_corrected(i);
    SMODS_cost = 2.2 + table.Pool_cost(i)*distance;

    p1 = table.p1(i)/100;
    p2 = table.p2(i)/100;
    p3 = table.p3(i)/100;
    p4 = table.p4(i)/100;
    p5 = table.p5(i)/100;
    p6 = table.p6(i)/100;

    if (table.Reference(i) == 'Public transit: ${q://QID56/ChoiceGroup/SelectedChoices}' && table.transit_mode(i) == 'Bus')

        % Scenario 1
        CE = walk_time * transit_walk * 1.1 + wait_time * 15 + transit_time * table.Bus_ref_1(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p1,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f1 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p1,R);

        % Scenario 2
        CE = walk_time * transit_walk * 1.1 + wait_time * 15 + transit_time * table.Bus_ref_2(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p2,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f2 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p2,R);
        
        % Scenario 3      
        CE = walk_time * transit_walk * 0.9 + wait_time * 1 + transit_time * table.Bus_ref_3(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 4 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 5 + wait_time * 6 + pool_time * 6.5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p3,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f3 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p3,R);
        
        % Scenario 4     
        CE = walk_time * transit_walk * 0.9 + wait_time * 1 + transit_time * table.Bus_ref_4(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 6 + wait_time * 5 + pool_time * 7 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 9 + wait_time * 4 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p4,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f4 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p4,R);
        
        % Scenario 5
        CE = walk_time * transit_walk + wait_time * 8 + transit_time * table.Bus_ref_5(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 6 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 2 + wait_time * 3 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p5,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f5 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p5,R);
        
        % Scenario 6
        CE = walk_time * transit_walk + wait_time * 5 + transit_time * table.Bus_ref_6(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 5 + wait_time * 6 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 2 + pool_time * 3.5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p6,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f6 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p6,R);

    elseif (table.Reference(i) == 'Public transit: ${q://QID56/ChoiceGroup/SelectedChoices}' && table.transit_mode(i) == 'Subway (T)')

        % Scenario 1
        CE = walk_time * transit_walk * 1.1 + wait_time * 9 + transit_time * table.Subway_ref_1(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p1,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f1 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p1,R);
        
        % Scenario 2
        CE = walk_time * transit_walk * 1.1 + wait_time * 9 + transit_time * table.Subway_ref_2(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p2,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f2 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p2,R);
        
        % Scenario 3      
        CE = walk_time * transit_walk * 0.9 + wait_time * 1 + transit_time * table.Subway_ref_3(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 4 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 5 + wait_time * 6 + pool_time * 8 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p3,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f3 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p3,R);
        
        % Scenario 4  
        CE = walk_time * transit_walk * 0.9 + wait_time * 1 + transit_time * table.Subway_ref_4(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 6 + wait_time * 5 + pool_time * 7 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 9 + wait_time * 4 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p4,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f4 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p4,R);
        
        % Scenario 5
        CE = walk_time * transit_walk + wait_time * 4 + transit_time * table.Subway_ref_5(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 6 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 2 + wait_time * 3 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p5,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f5 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p5,R);
        
        % Scenario 6
        CE = walk_time * transit_walk + wait_time * 5 + transit_time * table.Subway_ref_6(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 5 + wait_time * 6 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 2 + pool_time * 3.5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p6,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f6 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p6,R);

    elseif (table.Reference(i) == 'Public transit: ${q://QID56/ChoiceGroup/SelectedChoices}' && table.transit_mode(i) == 'Commuter rail')

        % Scenario 1
        CE = walk_time * transit_walk * 1.1 + wait_time * 15 + transit_time * table.Rail_ref_1(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 4 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 3 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p1,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f1 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p1,R);
        
        % Scenario 2
        CE = walk_time * transit_walk * 1.1 + wait_time * 20 + transit_time * table.Rail_ref_2(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 4 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 3 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p2,CE,i); 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f2 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p2,R);
        
        % Scenario 3      
        CE = walk_time * transit_walk * 0.9 + wait_time * 2 + transit_time * table.Rail_ref_3(i) + price * transit_cost + ASC_transit;        
        u1 = walk_time * 7 + wait_time * 4 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 5 + wait_time * 6 + pool_time * 6.5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p3,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f3 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p3,R);
        
        % Scenario 4      
        CE = walk_time * transit_walk * 0.9 + wait_time * 2 + transit_time * table.Rail_ref_4(i) + price * transit_cost + ASC_transit;        
        u1 = walk_time * 6 + wait_time * 5 + pool_time * 7 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 9 + wait_time * 4 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p4,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f4 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p4,R);
        
        % Scenario 5
        CE = walk_time * transit_walk + wait_time * 10 + transit_time * table.Rail_ref_5(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 6 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 2 + wait_time * 3 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p5,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f5 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p5,R);
        
        % Scenario 6
        CE = walk_time * transit_walk + wait_time * 10 + transit_time * table.Rail_ref_6(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 5 + wait_time * 6 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 2 + pool_time * 3.5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p6,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f6 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p6,R);

    elseif (table.Reference(i) == 'Exclusive ridesharing')

        exc_cost = 2.2 + table.Exc_cost(i) * table.Distance_corrected(i) + 0.42 * table.Exc_driver_wait2(i);

        % Scenario 1
        % Certainty equivalent - Sure prospect with certain alternative travel options
        CE = walk_time * 0 + wait_time * 9 + exc_time * table.Distance_corrected(i) * 5 + price * table.Exc_ref_1(i) + ASC_exclusive;
        u1 = walk_time * 5 + wait_time * 2 + pool_time * 3.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 4 + wait_time * 1 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p1,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f1 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p1,R);
        
        % Scenario 2
        CE = walk_time * 0 + wait_time * 9 + exc_time * table.Exc_ref_2(i) + price * exc_cost * 1.2 + ASC_exclusive;
        u1 = walk_time * 2 + wait_time * 4 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 3 + pool_time * 3.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p2,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f2 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p2,R);
        
        % Scenario 3
        CE = walk_time * 0 + wait_time * 1 + exc_time * table.Exc_ref_3(i) + price * exc_cost * 0.8 + ASC_exclusive;
        u1 = walk_time * 7 + wait_time * 4 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 5 + wait_time * 6 + pool_time * 6.5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p3,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f3 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p3,R);
        
        % Scenario 4
        CE = walk_time * 0 + wait_time * 1 + exc_time * table.Exc_ref_4(i) + price * exc_cost * 0.8 + ASC_exclusive;
        u1 = walk_time * 6 + wait_time * 5 + pool_time * 7 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 9 + wait_time * 4 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p4,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f4 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p4,R);
        
        % Scenario 5
        CE = walk_time * 0 + wait_time * 5 + exc_time * table.Exc_ref_5(i) + price * exc_cost + ASC_exclusive;
        u1 = walk_time * 7 + wait_time * 6 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 2 + wait_time * 3 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p5,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f5 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p5,R);
        
        % Scenario 6
        CE = walk_time * 0 + wait_time * 5 + exc_time * table.Exc_ref_6(i) + price * exc_cost + ASC_exclusive;
        u1 = walk_time * 5 + wait_time * 6 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 2 + pool_time * 3.5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p6,CE,i);
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f6 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p6,R);
              
    end
    
    x(5) = x(5)/100;
    F = [f1; f2; f3; f4; f5; f6; reg.*[sqrt(norm(x-lb)); sqrt(norm(x-ub))]];
    
end
