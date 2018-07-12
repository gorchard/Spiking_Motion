function output_events = implement_refraction_all(input_events, filters)
%function TDevts = ImplementRefractionAll(input_events, filters)
%   input_events    -> a struct of input events
%   filters         -> a 2D cell array of 2D filter delays

%find how many speeds and directions there are
[num_speeds, num_directions] = size(filters);
max_delay = zeros(1, num_speeds);

%find the maximum delay in the filters for each speed
for speed_index = 1:num_speeds
    max_delay(speed_index) = 0;
    for direction_index = 1:num_directions
        max_delay(speed_index) = max(max_delay(speed_index), max(filters{speed_index, direction_index}(:)));
    end
end


num_spikes = zeros(1, num_speeds); %find the number of output events for each speed
output_events = cell(1,num_speeds); %create a cell to hold output events for each speed
for speed_index = 1:num_speeds
    %implement a refraction period equal to the longest filter delay for each speed
    output_events{speed_index} = implement_refraction(input_events, max_delay(speed_index)); 
    
    %find the number of spikes in the refracted spike stream
    num_spikes(speed_index) = length(output_events{speed_index}.ts);
end

%find the longest of the refracted spike streams
max_num_spikes = max(num_spikes);

%make all the spikes streams the same size by padding with "Inf"
for speed_index = 1:num_speeds
    output_events{speed_index}.ts((end+1):max_num_spikes) = inf;
    output_events{speed_index}.x((end+1):max_num_spikes) = inf;
    output_events{speed_index}.y((end+1):max_num_spikes) = inf;
end