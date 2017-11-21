function output_spikes = spike_layer_2(input_spikes, filters, spike_threshold, decay_rate, reset_potential)
%function output_spikes = spike_layer_2(input_spikes, filters, spike_threshold, decay_rate, reset_potential)
%   input_spikes      -> The spikes to be filtered
%   filters           -> A 2D struct of size [num_speeds, num_directions] of 2D filters [x,y] to be used
%   spike_threshold   -> Integrate and Fire Threshold
%   decay_rate        -> Decay Rate in units per us
%   reset_potential   -> Potential immediately after reset
%
%The function relies on filters being spaced by 45 degrees in direction and by sqrt(2) in speed.

%% Constants
[num_speeds, num_directions] = size(filters); %how many speeds and directions are we using

filter_size = size(filters{1,1}); %how large is each filter (assume they're all the same size)
assert(rem(filter_size(1),2) == 1 & rem(filter_size(2),2) == 1); % check that the filter size is odd
filter_half_size = floor(filter_size./2);

dim_x = max(input_spikes.x(~isinf(input_spikes.x)));
dim_y = max(input_spikes.y(~isinf(input_spikes.y)));

synaptic_weight = 100/(filter_size(1)*filter_size(2)); %what is the weight of each synapse? (make the sum of the weights 100)
num_spikes_total = length(input_spikes.ts); %how many input spikes are there?
num_spikes_per_batch = 10000; %how many spikes to process per batch?
num_batches = ceil(num_spikes_total/num_spikes_per_batch); %how many batches must we process?

%% Organize incoming spikes 
%this struct will hold the chronologically (after synapse delays) ordered spikes
ordered_spikes.x = [];
ordered_spikes.y = [];
ordered_spikes.ts = [];
ordered_spikes.sp = [];
ordered_spikes.dir = [];

%this struct will hold the output spikes
output_spikes.x         = zeros(10000,1);
output_spikes.y         = zeros(10000,1);
output_spikes.ts        = zeros(10000,1);
output_spikes.sp        = zeros(10000,1);
output_spikes.dir       = zeros(10000,1);
%and the index where the next output spike should be placed
output_spike_number = 1; 

%arrays for remembering the last time a neuron was updated, and the neuron potential at that time
neuron_time = zeros(num_speeds, num_directions, dim_x, dim_y);
neuron_potential = zeros(num_speeds, num_directions, dim_x, dim_y);

%% Process
fprintf('Beginning processing of %i batches\n', num_batches);
tic %start a timer to keep track of progress
for batch_num = 1:num_batches %loop through the batches
    next_ordered_spike_index = length(ordered_spikes.x)+1; %find the index where the next ordered spike should be put
    num_spikes_this_batch = min(num_spikes_per_batch, length(input_spikes.x)); %find the batch size. Might be smaller than "num_spikes_per_batch" if too few spikes remain to be processed

    %to save time, pre-allocate memory for spikes which are about to be placed
    num_spikes_to_add = num_spikes_this_batch*3*filter_size(1)*filter_size(2); %how many spikes will be added to "ordered_spikes"?
    ordered_spikes.x    = [ordered_spikes.x; zeros(num_spikes_to_add, 1)];
    ordered_spikes.y    = [ordered_spikes.y; zeros(num_spikes_to_add, 1)];
    ordered_spikes.ts   = [ordered_spikes.ts; zeros(num_spikes_to_add, 1)];
    ordered_spikes.sp   = [ordered_spikes.sp; zeros(num_spikes_to_add, 1)];
    ordered_spikes.dir  = [ordered_spikes.dir; zeros(num_spikes_to_add, 1)];
    
    for x_filter = 1:filter_size(1)
        for y_filter = 1:filter_size(2)
            %add spikes from the current speed and direction
            ordered_spikes.x(next_ordered_spike_index:(next_ordered_spike_index+num_spikes_this_batch-1))  = input_spikes.x(1:num_spikes_this_batch) - x_filter + (filter_half_size(1) + 1);
            ordered_spikes.y(next_ordered_spike_index:(next_ordered_spike_index+num_spikes_this_batch-1))  = input_spikes.y(1:num_spikes_this_batch) - y_filter + (filter_half_size(2) + 1);
            for ordered_spike_index = 1:num_spikes_this_batch
                ordered_spikes.ts(next_ordered_spike_index+ordered_spike_index-1) = input_spikes.ts(ordered_spike_index) + filters{input_spikes.sp(ordered_spike_index),input_spikes.dir(ordered_spike_index)}(x_filter,y_filter);
            end
            ordered_spikes.sp(next_ordered_spike_index:(next_ordered_spike_index+num_spikes_this_batch-1)) = input_spikes.sp(1:num_spikes_this_batch);
            ordered_spikes.dir(next_ordered_spike_index:(next_ordered_spike_index+num_spikes_this_batch-1)) = input_spikes.dir(1:num_spikes_this_batch);
            next_ordered_spike_index = next_ordered_spike_index+num_spikes_this_batch;
            
            %add spikes to the neighbouring direction at a higher speed
            ordered_spikes.x(next_ordered_spike_index:(next_ordered_spike_index+num_spikes_this_batch-1))  = input_spikes.x(1:num_spikes_this_batch) - x_filter + (filter_half_size(1) + 1);
            ordered_spikes.y(next_ordered_spike_index:(next_ordered_spike_index+num_spikes_this_batch-1))  = input_spikes.y(1:num_spikes_this_batch) - y_filter + (filter_half_size(2) + 1);
            for ordered_spike_index = 1:num_spikes_this_batch
                ordered_spikes.ts(next_ordered_spike_index+ordered_spike_index-1) = input_spikes.ts(ordered_spike_index) + filters{input_spikes.sp(ordered_spike_index),input_spikes.dir(ordered_spike_index)}(x_filter,y_filter);
            end
            ordered_spikes.sp(next_ordered_spike_index:(next_ordered_spike_index+num_spikes_this_batch-1)) = input_spikes.sp(1:num_spikes_this_batch) + 1;
            ordered_spikes.dir(next_ordered_spike_index:(next_ordered_spike_index+num_spikes_this_batch-1)) = input_spikes.dir(1:num_spikes_this_batch) - 1;
            next_ordered_spike_index = next_ordered_spike_index+num_spikes_this_batch;
            
            %add spikes to the other neighbouring direction at a higher speed
            ordered_spikes.x(next_ordered_spike_index:(next_ordered_spike_index+num_spikes_this_batch-1))  = input_spikes.x(1:num_spikes_this_batch) - x_filter + (filter_half_size(1) + 1);
            ordered_spikes.y(next_ordered_spike_index:(next_ordered_spike_index+num_spikes_this_batch-1))  = input_spikes.y(1:num_spikes_this_batch) - y_filter + (filter_half_size(2) + 1);
            for ordered_spike_index = 1:num_spikes_this_batch
                ordered_spikes.ts(next_ordered_spike_index+ordered_spike_index-1) = input_spikes.ts(ordered_spike_index) + filters{input_spikes.sp(ordered_spike_index),input_spikes.dir(ordered_spike_index)}(x_filter,y_filter);
            end
            ordered_spikes.sp(next_ordered_spike_index:(next_ordered_spike_index+num_spikes_this_batch-1)) = input_spikes.sp(1:num_spikes_this_batch) + 1;
            ordered_spikes.dir(next_ordered_spike_index:(next_ordered_spike_index+num_spikes_this_batch-1)) = input_spikes.dir(1:num_spikes_this_batch) + 1;
            next_ordered_spike_index = next_ordered_spike_index+num_spikes_this_batch;
        end
    end
    
    % find the indices of spikes which are out of range (in speed or x,y location)
    valid_indices = (ordered_spikes.x > 0) & (ordered_spikes.y >0) & (ordered_spikes.x <= dim_x) & (ordered_spikes.y <= dim_y) & (ordered_spikes.sp <= num_speeds) & (ordered_spikes.sp > 0);
    ordered_spikes.ts        = ordered_spikes.ts(valid_indices);
    ordered_spikes.x         = ordered_spikes.x(valid_indices);
    ordered_spikes.y         = ordered_spikes.y(valid_indices);
    ordered_spikes.sp        = ordered_spikes.sp(valid_indices);
    ordered_spikes.dir       = ordered_spikes.dir(valid_indices);
    
    % Sort the spikes in chronological order
    [ordered_spikes.ts, order_index]            = sort(ordered_spikes.ts);
    ordered_spikes.x                            = ordered_spikes.x(order_index);
    ordered_spikes.y                            = ordered_spikes.y(order_index);
    ordered_spikes.sp                           = ordered_spikes.sp(order_index);
    ordered_spikes.dir                          = ordered_spikes.dir(order_index);
    
    %Take care of cyclic wrap around for directions
    ordered_spikes.dir(ordered_spikes.dir == 0) = num_directions;
    ordered_spikes.dir(ordered_spikes.dir == (num_directions+1)) = 1;
    
    %at this point we have a list of chronologically ordered spikes in  "ordered_spikes" which already take into account the synaptic delays
    
    %% Simulate the Neurons
    ordered_spike_index = 1;
    if ~isempty(ordered_spikes.ts) %check that there are spikes to process
        
        %we stop looping when the next input spike has a timestamp smaller than the current "ordered_spikes" (in which case we need to add more input spikes to ordered spikes before we continue processing)
        %or when we run out of ordered spikes to process
        while (ordered_spikes.ts(ordered_spike_index) < input_spikes.ts(num_spikes_this_batch)) && (ordered_spike_index < length(ordered_spikes.ts))
            
            %%perform the neuron update for the current ordered spike
            
            %implement decay towards 0
            if neuron_potential(ordered_spikes.sp(ordered_spike_index), ordered_spikes.dir(ordered_spike_index), ordered_spikes.x(ordered_spike_index), ordered_spikes.y(ordered_spike_index)) >= 0
                neuron_potential(ordered_spikes.sp(ordered_spike_index), ordered_spikes.dir(ordered_spike_index), ordered_spikes.x(ordered_spike_index), ordered_spikes.y(ordered_spike_index)) = max((neuron_potential(ordered_spikes.sp(ordered_spike_index), ordered_spikes.dir(ordered_spike_index), ordered_spikes.x(ordered_spike_index), ordered_spikes.y(ordered_spike_index)) - (ordered_spikes.ts(ordered_spike_index) - neuron_time(ordered_spikes.sp(ordered_spike_index), ordered_spikes.dir(ordered_spike_index), ordered_spikes.x(ordered_spike_index), ordered_spikes.y(ordered_spike_index)))*decay_rate(ordered_spikes.sp(ordered_spike_index))), 0);
            else
                neuron_potential(ordered_spikes.sp(ordered_spike_index), ordered_spikes.dir(ordered_spike_index), ordered_spikes.x(ordered_spike_index), ordered_spikes.y(ordered_spike_index)) = min((neuron_potential(ordered_spikes.sp(ordered_spike_index), ordered_spikes.dir(ordered_spike_index), ordered_spikes.x(ordered_spike_index), ordered_spikes.y(ordered_spike_index)) + (ordered_spikes.ts(ordered_spike_index) - neuron_time(ordered_spikes.sp(ordered_spike_index), ordered_spikes.dir(ordered_spike_index), ordered_spikes.x(ordered_spike_index), ordered_spikes.y(ordered_spike_index)))*decay_rate(ordered_spikes.sp(ordered_spike_index))), 0);
            end
            
            %add the synaptic weight
            neuron_potential(ordered_spikes.sp(ordered_spike_index), ordered_spikes.dir(ordered_spike_index), ordered_spikes.x(ordered_spike_index), ordered_spikes.y(ordered_spike_index)) = synaptic_weight + neuron_potential(ordered_spikes.sp(ordered_spike_index), ordered_spikes.dir(ordered_spike_index), ordered_spikes.x(ordered_spike_index), ordered_spikes.y(ordered_spike_index));
            
            %check if threshold has been exceeded
            if neuron_potential(ordered_spikes.sp(ordered_spike_index), ordered_spikes.dir(ordered_spike_index), ordered_spikes.x(ordered_spike_index), ordered_spikes.y(ordered_spike_index)) > spike_threshold
                
                %reset all neuron potentials at this location (mutual inhibition)
                neuron_potential(:, :, ordered_spikes.x(ordered_spike_index), ordered_spikes.y(ordered_spike_index)) = reset_potential; % Inhibit all other neurons at the same location
                
                %record the time at which we have updated these neurons
                neuron_time(:, :, ordered_spikes.x(ordered_spike_index), ordered_spikes.y(ordered_spike_index)) = ordered_spikes.ts(ordered_spike_index);
                
                %add the spike to our output spikes
                output_spikes.x(output_spike_number)    = ordered_spikes.x(ordered_spike_index);
                output_spikes.y(output_spike_number)    = ordered_spikes.y(ordered_spike_index);
                output_spikes.ts(output_spike_number)   = ordered_spikes.ts(ordered_spike_index);
                output_spikes.sp(output_spike_number)   = ordered_spikes.sp(ordered_spike_index);
                output_spikes.dir(output_spike_number)  = ordered_spikes.dir(ordered_spike_index);
                output_spike_number = output_spike_number + 1;
                
                %check if our output spike index exceeds the size of "output_spikes". 
                %if so, allocate space for another 10000 output spikes
                if output_spike_number > length(output_spikes.x)
                    output_spikes.x         = [output_spikes.x; zeros(10000,1)];
                    output_spikes.y         = [output_spikes.y; zeros(10000,1)];
                    output_spikes.ts        = [output_spikes.ts; zeros(10000,1)];
                    output_spikes.sp        = [output_spikes.sp; zeros(10000,1)];
                    output_spikes.dir       = [output_spikes.dir; zeros(10000,1)];
                end
            end
            
            %record the time at which we have updated the current neuron
            neuron_time(ordered_spikes.sp(ordered_spike_index), ordered_spikes.dir(ordered_spike_index), ordered_spikes.x(ordered_spike_index), ordered_spikes.y(ordered_spike_index)) = ordered_spikes.ts(ordered_spike_index);
            
            %increment the index to process the next "ordered_spike"
            ordered_spike_index = ordered_spike_index+1;
        end

        %remove all the spikes we have processed thus far from "ordered_spikes"
        ordered_spike_index = ordered_spike_index-1;
        ordered_spikes.ts(1:ordered_spike_index)    = [];
        ordered_spikes.x(1:ordered_spike_index)     = [];
        ordered_spikes.y(1:ordered_spike_index)     = [];
        ordered_spikes.dir(1:ordered_spike_index)   = [];
        ordered_spikes.sp(1:ordered_spike_index)    = [];
        
        %remove the spikes we've just used from "input_spikes"
        input_spikes.ts(1:num_spikes_this_batch)    = [];
        input_spikes.x(1:num_spikes_this_batch)     = [];
        input_spikes.y(1:num_spikes_this_batch)     = [];
        input_spikes.dir(1:num_spikes_this_batch)   = [];
        input_spikes.sp(1:num_spikes_this_batch)    = [];
    
    end
    fprintf('Layer 2 Loop %i finished at %f seconds\n', batch_num, toc);
    
    
end

if ~isempty(output_spikes)
    %remove any spikes for which space was allocated by not used
    output_spikes.x(output_spike_number:end)    = [];
    output_spikes.y(output_spike_number:end)    = [];
    output_spikes.ts(output_spike_number:end)   = [];
    output_spikes.sp(output_spike_number:end)   = [];
    output_spikes.dir(output_spike_number:end)  = [];
    
    %sort the output spikes by time (they should already be sorted, this shouldn't have any effect)
    [output_spikes.ts, order_index] = sort(output_spikes.ts);
    output_spikes.x                 = output_spikes.x(order_index);
    output_spikes.y                 = output_spikes.y(order_index);
    output_spikes.sp                = output_spikes.sp(order_index);
    output_spikes.dir               = output_spikes.dir(order_index);
    
    %give an arbitrary to the spikes. Mostly for compatibility with other functions
    output_spikes.p                 = ones(size(output_spikes.dir));
end
