%NewStream = combine_streams(Stream1, Stream2, newFieldName, newFieldVal1, newFieldVal2)
%   Combines two event streams (Stream1 and Stream2) into a single stream
%   (NewStream)
%   
%   For example, two TD events streams (TD1 and TD2) can be combined for
%   testing noise, where TD1 is signal and TD2 is noise
%
% TAKES IN:
% 'Stream1'
%   A struct of events. Typical format is:
%       Stream1.x =  pixel X locations
%       Stream1.y =  pixel Y locations
%       Stream1.p =  event polarity
%       Stream1.ts = event timestamps, typically in microseconds 
% 
% 'Stream2'
%   A struct of events with the same format as 'Stream1'
%   
% 'newFieldName'
%   An optional argument to add a new field to the output struct which
%   tracks which input struct each event is from
% 
% 'newFieldVal1'
%   An optional argument assigning a particular value to the
%   'newFieldName' field for events originating from 'Stream1'
% 
% 'newFieldVal2'
%   An optional argument assigning a particular value to the
%   'newFieldName' field for events originating from 'Stream2'
% 
% RETURNS:
% 'NewStream'
%   A struct with the same fields as 'Stream1', plus a new field
%   ('newFieldName') if the 'newFieldName' argument is used. The events in
%   'NewStream' will be sorted in chronological order (i.e. by the order of
%   the field 'ts').
% 
% 
% EXAMPLE USES:
% TDcombined = CombineStreams(TD1, TD2);
% % TDcombined now contains all events from TD1 and TD2.
% 
% 
% TDcombined = CombineStreams(TD1, TD2, 'source', 1, 2);
% % TDcombined now contains all events from TD1 and TD2
% % TDcombined.source = 1 if the event is from TD1
% % TDcombined.source = 2 if the event is from TD2
% 
% 
% written by Garrick Orchard - June 2014
% garrickorchard@gmail.com

function Stream1 = combine_streams(Stream1, Stream2, newFieldName, newFieldVal1, newFieldVal2)

% check the field names
fieldnames = fields(Stream1);
fieldnames2 = fields(Stream2);

% make sure they are all row vectors (not column vectors)
for i = 1:length(fieldnames)
    [numRows,~] = size(Stream1.(fieldnames{i}));
    if numRows>1
        Stream1.(fieldnames{i}) = Stream1.(fieldnames{i})';
    end
end

for i = 1:length(fieldnames2)
    [numRows,~] = size(Stream2.(fieldnames2{i}));
    if numRows>1
        Stream2.(fieldnames2{i}) = Stream2.(fieldnames2{i})';
    end
end



if exist('newFieldName', 'var')
    Stream1.(newFieldName)  =  [newFieldVal1*ones(size(Stream1.(fieldnames{1}))), newFieldVal2*ones(size(Stream2.(fieldnames{1})))];
end

if ~isempty(Stream2)
    for i = 1:length(fieldnames)
        if isfield(Stream2, fieldnames{i})
            Stream1.(fieldnames{i})  =  [Stream1.(fieldnames{i}), Stream2.(fieldnames{i})];
        else
            Stream1.(fieldnames{i})  =  [Stream1.(fieldnames{i}), zeros(size(Stream2.(fieldnames2{1})))];
        end
    end
end

if isfield(Stream1, 'TimeStamp') || isfield(Stream1, 'ts')
    Stream1 = sort_order(Stream1);
end