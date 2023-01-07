function fuel_mass_remain = FuelMassCalculation(fuel_mass,t,fuel_rate)
% Return the remained fuel mass of each node of the rocket at time t

%% find the starting and the end node of the second stage of the rocket 
ii = length(fuel_mass);
while(fuel_mass(ii)==0)
    ii = ii-1;
end
node_end_2nd_stage = ii;
while(fuel_mass(ii)>0)
    ii = ii-1;
end
node_start_2nd_stage = ii+1;


total_fuel_mass = sum(fuel_mass(node_start_2nd_stage:node_end_2nd_stage));
fuel_mass_remain = fuel_mass;
if total_fuel_mass < t*fuel_rate  % in this case, the fuel of the second stage has been exhausted
%     msgbox('The fuel of the second stage has been exhausted.');
    fuel_mass(node_start_2nd_stage:node_end_2nd_stage) = 0;
    fuel_mass_remain = fuel_mass;
else
    % in this case, the fuel of the second stage has not been exhausted
    N_node_2nd_stage = node_end_2nd_stage-node_start_2nd_stage+1; % total number of the nodes of the second stage
    temp = fuel_mass(node_start_2nd_stage:node_end_2nd_stage,1);
    taxis = zeros(N_node_2nd_stage,1);
    for ii = 1:N_node_2nd_stage
        taxis(ii,1) = sum(temp(1:ii))/fuel_rate;
    end
    temp = taxis-t;
    ind = find(temp>=0,1); % find the node that the fuel is being consumed at time t
    if ind>1
        fuel_mass_remain(node_start_2nd_stage:node_start_2nd_stage+ind-2,1) = 0;
        fuel_mass_remain(node_start_2nd_stage+ind-1,1) = fuel_mass(node_start_2nd_stage+ind-1,1)-fuel_rate*(t-taxis(ind-1));
    else
        fuel_mass_remain(node_start_2nd_stage,1) = fuel_mass(node_start_2nd_stage,1)-fuel_rate*t;
    end
end
end