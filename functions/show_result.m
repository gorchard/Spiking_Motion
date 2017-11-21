function vid = show_result(varargin)
% vid = ShowTD(TD, TPF, FL, time_span, Scale)
%   Shows a video of Temporal Difference (TD) events and returns a video
%   object which can be saved to AVI using the 'SaveMovie' function.  
%   All arguments except TD are optional.
%
% TAKES IN:
%   'TD'
%       A struct of TD events with format:
%           TD.x =  pixel X locations
%           TD.y =  pixel Y locations
%           TD.p =  event polarity
%           TD.ts = event timestamps in microseconds
% 
%   'TPF'
%       Time Per Frame (TPF) is an optional argument specifying the
%       time-spacing between the start of subsequent frames (basically the
%       frame rate for the video). Defaults to 24FPS, which is a Time Per
%       Frame (TPF) of 1/24 seconds. 
% 
%   'FL'
%       Frame Length (FL) is an optional arguments specifying the time-span
%       of data to show per frame as a fraction of TPF. Defaults to TPF seconds. If FL<1,
%       then not all data in a sequence will be shown. If FL>1 then some
%       data will be repeated in subsequent frames.
% 
%   'time_span' = [Tstart,Tstop]
%       An optional argument specifying at which point in time the playback
%       should start (Tstart) and stop (Tstop). If time_span is not
%       specified, the entire recording will be shown by default.
% 
% 
% RETURNS:
%    'vid' 
%       A video object which can be saved to AVI using the 'SaveMovie'
%       function.
%
% 
% written by Garrick Orchard - June 2014
% garrickorchard@gmail.com

close all
s = warning('off','images:imshow:magnificationMustBeFitForDockedFigure');
timeconst = 1e-6;
TD = varargin{1};

%TD.p = round(TD.p - min(TD.p) + 1);

%FPS is 1/TPF
if nargin > 1
    if isempty(varargin{2})
        FPS = 24;
    else
        FPS = 1/varargin{2};
    end
else
    FPS = 24;
end

%FL is overlap
if nargin >2
    if isempty(varargin{3})
        Overlap = 1;
    else
        Overlap = varargin{3};
    end
else
    Overlap = 1;
end

if nargin > 3
    if isempty(varargin{4})
        Tmin = 1;
        Tmax = length(TD.ts);
    else
        if(varargin{4}(1) == -1)
            Tmin = 1;
        else
            Tmin = find(TD.ts>varargin{4}(1),1);
        end
        if(varargin{4}(2) == -1)
            Tmax = length(TD.ts);
        else
            Tmax = find(TD.ts>varargin{4}(2),1);
        end
        if isempty(Tmax)
            Tmax = length(TD.ts);
        end
    end
else
    Tmin = 1;
    Tmax = length(TD.ts);
end

FrameLength = 1/(FPS*timeconst);
t1 = TD.ts(Tmin) + FrameLength;
t2 = TD.ts(Tmin) + FrameLength*Overlap;

dim_x = max(TD.x);
dim_y = max(TD.y);
ImageBack = zeros(dim_y,dim_x*2,3);

axis image
i = Tmin;
cc.sp = jet(double(max(TD.sp))); %speed colors (jet for a scale)
cc.dir = hsv(double(max(TD.dir))); %direction colors (hsv because it wraps around)
Image = ImageBack;
k=1;
nFrames = ceil((TD.ts(Tmax)-TD.ts(Tmin))/FrameLength);
vid(1:nFrames) = struct('cdata', ImageBack, 'colormap', []);
while (i<Tmax)
    j=i;
    while ((TD.ts(j) < t2) && (j<Tmax))
        Image(TD.y(j), TD.x(j), :) = cc.dir(TD.dir(j),:);
        Image(TD.y(j), TD.x(j)+dim_x, :) = cc.sp(TD.sp(j),:);
        j = j+1;
    end
    while ((TD.ts(i) < t1) && (i<Tmax))
        i = i+1;
    end
    imshow(Image, 'InitialMagnification', 'fit');
    text(dim_x/2, dim_y-10, 'Direction', 'HorizontalAlignment', 'center', 'color' , [1, 1, 1])
    text(dim_x/2+dim_x, dim_y-10, 'Speed', 'HorizontalAlignment', 'center', 'color' , [1, 1, 1])
    title(TD.ts(i))
    axis off
    drawnow();
%     input(' ')
    
    t2 = t1 + FrameLength*Overlap;
    t1 = t1 + FrameLength;
    vid(k).cdata = Image;
    Image = ImageBack;
    k=k+1;
end