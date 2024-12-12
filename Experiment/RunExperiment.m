%RunExperiment(sub,run,mode)

function RunExperiment(par,run,mode)

% % % for Kevin's testing
% % nargin = 3;
% % par = 0;
% % run = 2;
% % mode = 'noprojector';

%% Uncomment this if the screen won't start even after "clear all" and "close all" (image presentation timing will not be good)
%Screen('Preference', 'SkipSyncTests', 1);

%% Check PsychToolbox
AssertOpenGL();

%% Inputs
if nargin < 2
    error(sprintf('%s\nToo few inputs!',help(mfilename)))
end
if ~exist('mode', 'var')
    mode = 'testing';
end
mode = lower(mode);
switch mode
    case 'testing'
        p.DEBUG = false;
        p.USE_PROJECTOR = true;
    case 'debug'
        p.DEBUG = true;
        p.USE_PROJECTOR = true;
    case 'noprojector'
        p.DEBUG = false;
        p.USE_PROJECTOR = false;
    case 'debug_noprojector'
        p.DEBUG = true;
        p.USE_PROJECTOR = false;
    otherwise
        error('Unknown mode!');
end

%% Parameters

%TR in seconds
p.TR = 1;

%paths
p.PATH.ORDERS_FOLDER = [pwd filesep 'Orders' filesep];
p.PATH.DATA_FOLDER = [pwd filesep 'Data'  filesep];
p.PATH.IMAGE = [pwd filesep 'Images' filesep];
p.PATH.ORDER = [p.PATH.ORDERS_FOLDER sprintf('PAR%02d_RUN%02d.xls*',par,run)];
p.PATH.DATA = [p.PATH.DATA_FOLDER sprintf('PAR%02d_RUN%02d_%s_%s',par,run,mode,get_timestamp)];

%trigger checking
p.TRIGGER.TIME_BEFORE_TRIGGER_MUST_START_LOOKING_SEC = 0.010; %should be less than TR
p.TRIGGER.TIME_BEFORE_TRIGGER_CAN_START_LOOKING_SEC = 0.500;
p.TRIGGER.TIME_AFTER_MISSED_TRIGGER_STOP_LOOKING_SEC = 0.005;

%misc
p.KEY.STOP = 27; %ESC key
p.KEY.TRIGGER = [53 84]; %5 and/or T
p.KEY.BUTTON_BOX = [49 50 51 52 82 71 66 89 97 98 99 100]; %1-4 top of key board, rgby, 1-4 numpad

%screen
p.SCREEN.NUMBER = 0;
p.SCREEN.RECT = []; %[] is full screen, [0 0 1920 1080] is the intended size as a window
p.SCREEN.EXPECTED_SIZE = [1080 1920]; %[height width]

%colours
p.SCREEN.BACKGROUND_COLOUR = [0 0 0];
p.SCREEN.TEXT_COLOUR = [255 255 255];

%stereo
p.SCREEN.STEREO_MODE = 1; %1 or 11 for shutter glasses (1 seems to work better)
p.SCREEN.PIPELINE_MODE = kPsychNeedFastBackingStore;
p.SCREEN.BUFFER_ID.LEFT = 0; %flip these if L/R is reversed
p.SCREEN.BUFFER_ID.RIGHT = 1;

%image
p.IMAGES.EXPECTED_WIDTH = 1920;
p.IMAGES.EXPECTED_HEIGHT = 1080;
p.IMAGES.VERTICAL_SHIFT = 0; %moves everything vertically, in units of PIXELS, positive is down and negative is up
p.IMAGES.RESIZE_FACTOR = 1; %loaded images are scaled by this value, also adjusts EXPECTED_WIDTH and EXPECTED_HEIGHT, 1.1=10% larger

% % % for Kevin's testing
% % Screen('Preference', 'SkipSyncTests', 1);
% % p.SCREEN.NUMBER = 1;
% % p.SCREEN.RECT = [0 0 1920 1080];

%% Prepare

%adjust expected dimensions
p.IMAGES.EXPECTED_WIDTH = round(p.IMAGES.EXPECTED_WIDTH * p.IMAGES.RESIZE_FACTOR);
p.IMAGES.EXPECTED_HEIGHT = round(p.IMAGES.EXPECTED_HEIGHT * p.IMAGES.RESIZE_FACTOR);

%make future calls faster
GetSecs;
KbCheck;

%create data folder if needed
if ~exist(p.PATH.DATA_FOLDER), mkdir(p.PATH.DATA_FOLDER);, end

%store git repo info
if exist('IsGitRepo','file')
      if ~IsGitRepo
         warning('This project does not appear to be part of a git repository. No git data will be saved.');
     elseif exist('GetGitInfo','file')
         d.GitInfo = GetGitInfo;
      end
 else
     warning('The "CulhamLab/Git-Version" repo has not been configured. Information about this project''s current repository status (version, etc.) will NOT be saved to the data file.');
 end

%load order
list = dir(p.PATH.ORDER);
if isempty(list)
    error('Could not locate order file: %s', p.PATH.ORDER)
elseif length(list)>1
    error('Multiple matches for order file: %s', p.PATH.ORDER);
else
    p.PATH.ORDER = [list.folder filesep list.name];
end
[~,~,d.order] = xlsread(p.PATH.ORDER);
d.order_headers = d.order(1,:);

%remove empty rows at end
row_last = find(cellfun(@(x) length(x)==1 && isnan(x), d.order(:,2)) == 0, 1, 'last');
d.order = d.order(1:row_last,:);

%create baseline frame
imgBaseline = uint8(zeros([p.IMAGES.EXPECTED_HEIGHT p.IMAGES.EXPECTED_WIDTH 3]));
imgBaseline(:,:,1) = p.SCREEN.BACKGROUND_COLOUR(1);
imgBaseline(:,:,2) = p.SCREEN.BACKGROUND_COLOUR(2);
imgBaseline(:,:,3) = p.SCREEN.BACKGROUND_COLOUR(3);
s.imgBaselineLeft = imgBaseline;
s.imgBaselineRight = imgBaseline;

%load images
[d,images] = load_images(d,p,s);

%% Volume Schedule

%COL 1: trial
%COL 2: image number to display, left view (0=no image)
%COL 3: image display duration in seconds (0=don't stop)
%COL 4: image number to prepare, left view (0=don't load)
%COL 5: order row
%COL 6: image number to display, right view (0=no image)
%COL 7: image number to prepare, right view (0=don't load)

%init schedule
d.sched = [];

for row = 2:size(d.order,1)
    if (d.order{row,3} / p.TR) ~= round(d.order{row,3} / p.TR)
        error('Events cannot be fraction of TR.')
    end
    
    if isempty(d.order{row,4}) && isempty(d.order{row,5})
        %NULL
        for v = 1:p.TR:d.order{row,3}
            d.sched(end+1,:) = [0 0 0 0 row 0];
        end
    else
        %not NULL
        filename_left = d.order{row,4};
        image_number_left = find(strcmp(d.image_filename_lookup,filename_left));
        if isempty(image_number_left)
            image_number_left = 0;
        end
        filename_right = d.order{row,5};
        image_number_right = find(strcmp(d.image_filename_lookup,filename_right));
        if isempty(image_number_right)
            image_number_right = 0;
        end
        % if any(strcmpi(p.MISC.CONDITIONS_WITH_UNLIMITED_DISPLAY_DURATION, d.order{row,2}))
            dur = inf;
        % else
        %     dur = p.DURATION.IMAGE_PRESENTATION_SECONDS;
        % end
        for v = 1:p.TR:d.order{row,3}
            image_number_left_this = image_number_left;
            image_number_right_this = image_number_right;
            if dur > p.TR
                dur_this = 0;
                dur = dur - p.TR;
            elseif dur <= 0
                dur_this = 0;
                image_number_left_this = 0;
                image_number_right_this = 0;
            else
                dur_this = dur;
                dur = 0;
            end
            d.sched(end+1,:) = [d.order{row,1} image_number_left_this dur_this 0 row image_number_right_this];
        end
    end
end

%fill COL 4/7
for row = 1:(size(d.sched,1)-1)
    if d.sched(row+1,2)>0
        d.sched(row,4) = d.sched(row+1,2);
    end
    if d.sched(row+1,6)>0
        d.sched(row,7) = d.sched(row+1,6);
    end
end

d.number_volumes = size(d.sched,1);

%% Try...
try

% Enable PROPixx RB3D Sequencer
if p.USE_PROJECTOR
    Datapixx('Open'); 
    Datapixx('EnableVideoStereoBlueline');
end

%% Screen

[s.win, s.rect] = Screen('OpenWindow', p.SCREEN.NUMBER, p.SCREEN.BACKGROUND_COLOUR, p.SCREEN.RECT, [], [], p.SCREEN.STEREO_MODE, [], p.SCREEN.PIPELINE_MODE);
if s.rect(1)~=0 || s.rect(2)~=0 || s.rect(3)~=p.SCREEN.EXPECTED_SIZE(2) || s.rect(4)~=p.SCREEN.EXPECTED_SIZE(1)
    error('Unexpected screen size! [%s]',num2str(s.rect))
end
HideCursor;

%% Clut

%set GPU CLUTs to linear
if p.USE_PROJECTOR
    Screen('LoadNormalizedGammaTable', s.win, linspace(0,1,256)'*[1,1,1]);
end

%% Make Image Textures

textures.baseline.left = Screen('MakeTexture', s.win, s.imgBaselineLeft);
textures.baseline.right = Screen('MakeTexture', s.win, s.imgBaselineRight);

for i = 1:length(images)
    textures.images(i).left = Screen('MakeTexture', s.win, images(i).left);
    textures.images(i).right = Screen('MakeTexture', s.win, images(i).right);
end

%% Wait for first trigger

message = sprintf('\nWaiting For Keys To Release',d.number_volumes,p.TR);
StereoDraw(s.win, p, textures, 0, 0, message);
Screen('DrawingFinished',s.win);
Screen('Flip',s.win);

%make sure nothing is pressed currently
while KbCheck, end

%first frame

message = sprintf('\nWaiting For Trigger (%d vol, TR %d)',d.number_volumes,p.TR);
StereoDraw(s.win, p, textures, 0, 0, message);
Screen('DrawingFinished',s.win);
Screen('Flip',s.win);

ready_to_draw = false;

%wait for first trigger
while 1
    %%%%KbWait; %more efficient way to wait for a key
    [keyIsDown, ~, keyCode] = KbCheck(-1); %get key(s)
    if keyIsDown
        if any(keyCode(p.KEY.TRIGGER))
            t0 = GetSecs;
            break
        elseif any(keyCode(p.KEY.STOP))
            error('Stop key was pressed.')
        end
    end
end

%% Run Volumes
vol_first = 1;
for v = vol_first:d.number_volumes
    %volume start time actual
    if v==vol_first
        d.volume_data(v).time_startActual = 0;
    else
        d.volume_data(v).time_startActual = GetSecs-t0;
    end
    
    %volume start time
    if v==vol_first | d.volume_data(v-1).recievedTrigger %is first vol OR prior vol recieved trigger
        d.volume_data(v).time_start = d.volume_data(v).time_startActual; %use actual time
    else %missed a trigger
        d.volume_data(v).time_start = d.volume_data(v-1).time_start + p.TR; %use expected trigger time
    end
    
    %start message
    fprintf('\nStarting volume %d/%d at %fsec (actual %fsec):\n',v,d.number_volumes,d.volume_data(v).time_start,d.volume_data(v).time_startActual);
    
    %info
    d.volume_data(v).trial = d.sched(v,1);
    d.volume_data(v).image_number_left = d.sched(v,2);
    d.volume_data(v).image_number_right = d.sched(v,6);
    if d.volume_data(v).image_number_left
        d.volume_data(v).image_filename_left = images(d.volume_data(v).image_number_left).filename;
    else
        d.volume_data(v).image_filename_left = nan;
    end
    if d.volume_data(v).image_number_right
        d.volume_data(v).image_filename_right = images(d.volume_data(v).image_number_right).filename;
    else
        d.volume_data(v).image_filename_right = nan;
    end
    if d.sched(v,3)
        d.volume_data(v).time_to_remove_image = d.sched(v,3);
    else
        d.volume_data(v).time_to_remove_image = nan;
    end
    d.volume_data(v).next_image_left = d.sched(v,4);
    d.volume_data(v).next_image_right = d.sched(v,7);
    d.volume_data(v).button_press = false;
    d.volume_data(v).button_press_time = nan;
    
    %add excel info
    row = d.sched(v,5);
    fprintf('-Condition:\n');
    for f = 1:length(d.order_headers)
        eval(sprintf('d.volume_data(v).condition.%s = d.order{row,f};',d.order_headers{f}));
        val = d.volume_data(v).condition.(d.order_headers{f});
        if isnumeric(val) || islogical(val)
            val = num2str(val);
        end
        fprintf('\t%s: %s\n', d.order_headers{f}, val);
    end
    
    %first draw
    if ~ready_to_draw
        if p.DEBUG
            message = sprintf('%d',v);
        else
            message = [];
        end
        StereoDraw(s.win, p, textures, d.volume_data(v).image_number_left, d.volume_data(v).image_number_right, message);
        Screen('DrawingFinished',s.win);
    end
    Screen('Flip',s.win);
    d.volume_data(v).time_first_draw = (GetSecs-t0) - d.volume_data(v).time_start;
    d.volume_data(v).time_second_draw = [];
    ready_to_draw = false;
    
    %will there be another draw?
    if isnan(d.volume_data(v).time_to_remove_image)
        done_draw_this_vol = true;
    else
        done_draw_this_vol = false;
    end
    
    %run volume events until it is time to check for trigger
    saved = false;
    nothing_more_to_do = false;
    while 1
        timeInVol = (GetSecs-t0) - d.volume_data(v).time_start;
        if (p.TR-timeInVol)<=p.TRIGGER.TIME_BEFORE_TRIGGER_MUST_START_LOOKING_SEC | (nothing_more_to_do && (p.TR-timeInVol)<=p.TRIGGER.TIME_BEFORE_TRIGGER_CAN_START_LOOKING_SEC)
            break
        elseif ~saved && nothing_more_to_do && (p.TR-timeInVol)>p.TRIGGER.TIME_BEFORE_TRIGGER_CAN_START_LOOKING_SEC
            %save if there is lots of time left and nothing else to do
            save(p.PATH.DATA,'d','p','s')
            saved=true;
        end
        
        if ~ready_to_draw
            if ~done_draw_this_vol
                %prep the next draw of this volume
                if p.DEBUG
                    message = sprintf('%d',v);
                else
                    message = [];
                end
                StereoDraw(s.win, p, textures, 0, 0, message);
                Screen('DrawingFinished',s.win);
                ready_to_draw = true;
            else
                %prep the next volume draw
                if p.DEBUG
                    message = sprintf('%d',v+1);
                else
                    message = [];
                end
                StereoDraw(s.win, p, textures, d.volume_data(v).next_image_left, d.volume_data(v).next_image_right, message);
                Screen('DrawingFinished',s.win);
                ready_to_draw = true;
            end
        end
        
        %nothing more to do?
        nothing_more_to_do = (ready_to_draw && done_draw_this_vol);
        
        %mid-volume draw
        if ~done_draw_this_vol && timeInVol>d.volume_data(v).time_to_remove_image
            if ~ready_to_draw
                if p.DEBUG
                    message = sprintf('%d',v+1);
                else
                    message = [];
                end
                StereoDraw(s.win, p, textures, 0, 0, message);
                Screen('DrawingFinished',s.win);
            end
            Screen('Flip',s.win);
            d.volume_data(v).time_second_draw = timeInVol;
            ready_to_draw = false;
            done_draw_this_vol = true;
        end
        
        %keys/button box
        [keyIsDown, ~, keyCode] = KbCheck(-1); %get key(s)
        if keyIsDown
            if any(keyCode(p.KEY.STOP))
                error('Stop key was pressed.') 
            elseif ~d.volume_data(v).button_press & any(keyCode(p.KEY.BUTTON_BOX))
                d.volume_data(v).button_press = true;
                d.volume_data(v).button_press_time = timeInVol;
                fprintf('-Button Box\n');
            end
        end
        
    end
    
    %look for trigger until time runs out
    d.volume_data(v).recievedTrigger = false; %default
    while 1
        timeInVol = (GetSecs-t0) - d.volume_data(v).time_start;
        
        if timeInVol>(p.TR+p.TRIGGER.TIME_AFTER_MISSED_TRIGGER_STOP_LOOKING_SEC)
            warning('No trigger was recieved. Continuing with expected timing...')
            break
        end
        
        [keyIsDown, ~, keyCode] = KbCheck(-1); %get key(s)
        if keyIsDown
            if any(keyCode(p.KEY.TRIGGER))
                d.volume_data(v).recievedTrigger = true;
                fprintf('~~~~~~~~~~~~~~~~~TRIGGER RECIEVED~~~~~~~~~~~~~~~~~\n')
                break
            elseif any(keyCode(p.KEY.STOP))
                error('Stop key was pressed.') 
            elseif ~d.volume_data(v).button_press & any(keyCode(p.KEY.BUTTON_BOX))
                d.volume_data(v).button_press = true;
                d.volume_data(v).button_press_time = timeInVol;
                fprintf('-Button Box\n');
            end
        end
    end
    
    %end of volume timing
    d.volume_data(v).time_endActual = GetSecs-t0;
    d.volume_data(v).volDuration = d.volume_data(v).time_endActual - d.volume_data(v).time_start;
    d.volume_data(v).volDurationActual = d.volume_data(v).time_endActual - d.volume_data(v).time_startActual;
    fprintf('-duration: %f seconds\n',d.volume_data(v).volDuration)
    
end

%% Report button press
d.number_button_presses = sum([d.volume_data.button_press]);
fprintf('\n\nButton Presses: %d\n\n%s\n\n',d.number_button_presses,num2str(find([d.volume_data.button_press])));

%% Cleanup
sca
sca
ShowCursor;
save(p.PATH.DATA,'d','p','s')

%PROPixx complete
if p.USE_PROJECTOR
    % Set the PROPixx back to normal sequencer
    Datapixx('SetPropixxDlpSequenceProgram', 0);
    Datapixx('RegWrRd');

    % Close PROPixx connection
    Datapixx('Close');
end


%% Done
disp Done!



%% Catch
catch err
    err
    sca
    sca
    ShowCursor;
    clear images
    save
    
    %PROPixx complete
    if p.USE_PROJECTOR
        % Set the PROPixx back to normal sequencer
        Datapixx('SetPropixxDlpSequenceProgram', 0);
        Datapixx('RegWrRd');
        
        % Close PROPixx connection
        Datapixx('Close');
    end
    
    rethrow(err)
end

function [timestamp] = get_timestamp
c = round(clock);
timestamp = sprintf('%d-%d-%d_%d-%d_%d',c([4 5 6 3 2 1]));

function [d,images] = load_images(d,p,s)

%get all image names
img_names = d.order(2:end,4:5);
img_names = img_names(:);

%remove non-string names (nan etc)
ind_valid = find(cellfun(@isstr,img_names));
img_names = img_names(ind_valid);

%reduce to unique names
d.image_filename_lookup = unique(img_names);

%convert to cell format (might not be cell if there is only one unique image)
if ~iscell(d.image_filename_lookup)
    d.image_filename_lookup = {d.image_filename_lookup};
end

%load images
for f = 1:length(d.image_filename_lookup)
    images(f).filename = d.image_filename_lookup{f};
    image = imread([p.PATH.IMAGE images(f).filename]);
   
    %resize
    if p.IMAGES.RESIZE_FACTOR ~= 1
        image = imresize(image, p.IMAGES.RESIZE_FACTOR);
    end
    
    %trim image where larger than expected
    width = size(image,2);
    trim_width = width - p.IMAGES.EXPECTED_WIDTH;
    if trim_width > 0
        select = (1 + ceil(trim_width/2)) : (width - floor(trim_width/2));
        image = image(:,select,:);
    end
    height = size(image,1);
    trim_height = height - p.IMAGES.EXPECTED_HEIGHT;
    if trim_height > 0
        select = (1 + ceil(trim_height/2)) : (height - floor(trim_height/2));
        image = image(select,:,:);
    end
    
    %vertical shift
    if p.IMAGES.VERTICAL_SHIFT > 0
        sz = size(image);
        image = [zeros(p.IMAGES.VERTICAL_SHIFT,sz(2),3,'uint8'); image];
        max_y = min(sz(1) + p.IMAGES.VERTICAL_SHIFT, p.IMAGES.EXPECTED_HEIGHT);
        image = image(1:max_y, :, :);
    elseif p.IMAGES.VERTICAL_SHIFT < 0
        pix = abs(p.IMAGES.VERTICAL_SHIFT);
        sz = size(image);
        image = [image; zeros(pix,sz(2),3,'uint8')];
        max_y = min(sz(1) + pix, p.IMAGES.EXPECTED_HEIGHT);
        image = image(end-max_y+1:end, :, :);
    end
    
    %fill edges where larger than expected
    width = size(image,2);
    width_fill = p.IMAGES.EXPECTED_WIDTH - width;
    if width_fill > 0
        height = size(image,1);
        image = [repmat(reshape(p.SCREEN.BACKGROUND_COLOUR,[1 1 3]),[height ceil(width_fill/2)]) image repmat(reshape(p.SCREEN.BACKGROUND_COLOUR,[1 1 3]),[height floor(width_fill/2)])];
    end
    height = size(image,1);
    height_fill = p.IMAGES.EXPECTED_HEIGHT - height;
    if height_fill > 0
        width = size(image,2);
        image = [repmat(reshape(p.SCREEN.BACKGROUND_COLOUR,[1 1 3]),[ceil(height_fill/2) width]); image; repmat(reshape(p.SCREEN.BACKGROUND_COLOUR,[1 1 3]),[floor(height_fill/2) width])];
    end
    
    %check size
    [height,width,~] = size(image);
    if p.IMAGES.EXPECTED_WIDTH ~= width || p.IMAGES.EXPECTED_HEIGHT ~= height
        p.IMAGES
        height
        width
        error(sprintf('%s is not the correct size!',images(f).filename));
    end
    
    %convert greyscale format to RGB
    if size(image,3) == 1
        image = repmat(image,[1 1 3]);
    end
    
    %degamma'ing the images so they don't look as pale
        image = imadjust(image, [], [], 1/log(2.2));
        
    %copies for left and right view
    images(f).left = image;
    images(f).right = image;
    
end

function StereoDraw(window, param, textures, index_left, index_right, message)
has_msg = ~isempty(message) && ischar(message);

if index_left == 0
    left = textures.baseline.left;
else
    left = textures.images(index_left).left;
end

if index_right == 0
    right = textures.baseline.right;
else
    right = textures.images(index_right).right;
end
%set GPU CLUTs to linear
Screen('LoadNormalizedGammaTable', window, linspace(0,1,256)'*[1,1,1]);

Screen('SelectStereoDrawBuffer', window, param.SCREEN.BUFFER_ID.LEFT);
Screen('DrawTexture', window, left);
if has_msg
    DrawFormattedText(window, message, 0, 20, param.SCREEN.TEXT_COLOUR);
end
    
Screen('SelectStereoDrawBuffer', window, param.SCREEN.BUFFER_ID.RIGHT);
Screen('DrawTexture', window, right);
if has_msg
    DrawFormattedText(window, message, 0, 20, param.SCREEN.TEXT_COLOUR);
end