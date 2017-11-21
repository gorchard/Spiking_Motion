function [ filters ] = generate_filters( varargin )
%generate_filters(speed_vector, direction_vector, filter_size)
%   speed_vector        -> a vector of the speeds (pix/ms)
%   direction_vector    -> a vector of the directions (degrees)
%   filter_size         -> size of each filter (2 element vector)

speed_vector = varargin{1};
direction_vector = varargin{2};
filter_size = varargin{3};

%allocate a cell to hold the filters
filters = cell(length(speed_vector), length(direction_vector));

for speed_index = 1:length(speed_vector)
    for direction_index = 1:length(direction_vector)
        
        %allocate the filter size
        filters{speed_index, direction_index} = zeros(filter_size(1), filter_size(2));
        
        %compute each delay
        for x = 1:filter_size(1)
            for y = 1:filter_size(2)
                filters{speed_index, direction_index}(x,y) = -(x*cosd(direction_vector(direction_index)) + y*sind(direction_vector(direction_index)))./(speed_vector(speed_index)/(1e3));
            end
        end
        
        %subtract the smallest delay to keep delays as small as possible
        filters{speed_index, direction_index} = filters{speed_index, direction_index} - min(filters{speed_index, direction_index}(:));
        
    end
end