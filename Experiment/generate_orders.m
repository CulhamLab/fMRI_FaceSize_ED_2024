% Note that these orders are balanced across runs, but not across participants
% There is no need to use a specific number/increment of participant orders.

% Rerunning this script for the same participant(s) will result in indentical
% orders. Comment out "rng(par)" to instead make pseudorandom orders.

% Overview (for more detail, search for "RULE")
% -6 runs
% -each run contains 12 blocks: 3 dist x 4 rep
% -the 4 reps are divided into 2 in each half of the run
% -distances cannot follow themselves
% -in the 3x3 distance order matrix (no diag), 1 cell will be "1" in each run (fully balanced across 6 runs)
% -2 runs will start with each of the 3 distances
% -distance blocks include 6 trials (2 of each face size), 1s duration, 3s ITI
% -independently for each distance... the sizes follow themself exactly 2
%   times per run and follow each other size 2-3 times per run (x2 in 4 runs, x3 in 2 runs).
% -20s baselines at the start and end
% -20s vergence change period before each block (including the first)
% -There are 18 faces (enough for 3 blocks). These will occur once in each 1/4 of a run.
% -For each of the 9 conditions, faces will occur 0 or 1 times per run, and
%   between 1-3 times total (suboptimal balance due to limitations from
%   combined rules, 2-3 total would be optimal)
% -Each face will follow each other face 0-1 times per run and 1-2 times
%   total. Faces never follow themself.

% IMPORTANT
% This script has been written for the numbers above. Any changes to the 
% numbers may break counterbalancing, cause errors, or even cause the 
% script to run forever if there is no longer a solution.
%

%% 
clear all


%% Which orders to make?

participant_orders_to_make = 1:30;


%% Stims

stims.faces = arrayfun(@(x) sprintf("Face%02d", x), 1:18);
stims.faces_count = length(stims.faces);

stims.sizes = arrayfun(@(x) sprintf("S%02d", x), [20 25 30]);
stims.sizes_count = length(stims.sizes);

stims.distances = arrayfun(@(x) sprintf("D%02d", x), [51 64 80]);
stims.distances_count = length(stims.distances);

stims.baseline_distance = "D64";

stims.image_prefix_baseline = "Fixation";
stims.image_prefix_vergence_change = "Fixation";
stims.image_prefix_ITI = "Fixation";

stims.image_suffix_left_eye = "l";
stims.image_suffix_right_eye = "r";
stims.image_filetype = ".png";


%% Parameters

fol = ".\Orders\";

% durations (seconds)
durations.baseline_initial = 20;
durations.baseline_vergence_change = 20;
durations.stim = 1;
durations.ITI = 3;
durations.baseline_final = 20;


%% Prep

% adjust slashes for OS
fol = fol.replace("/",filesep);
fol = fol.replace("\",filesep);

% output folder should end with a slash
if ~fol.endsWith(filesep)
    fol = fol + filesep;
end

% make output folder?
if ~exist(fol, "dir")
    mkdir(fol)
end


%% Participant Loop...
for par = participant_orders_to_make
fprintf("Generating orders for PAR%02d...\n", par);

% check if RUN01 order already exists
fp = fol + sprintf("PAR%02d_RUN01.xlsx", par);
if exist(fp, "file")
    fprintf("\tRUN01 already exist. Skipping this participant!\n")
    continue
end

% use fixed randomization based on participant number so that results are
% consistent / reproducable
rng(par)


%% Part 1: Order of Distance Blocks Over the 6 Runs
% This process will always find a solution without needing to restart

fprintf("\tGenerating distance block orders...\n");

% RULE: 2 of the 6 runs will start with each of the 3 distances.
% To spread these out, there will be 1 in each half with no repeat
% between the 3rd and 4th runs.
%
% (randomly 1-3) (randomly 1-3)
% but never [. . X X . .]
%
run_first_distances = [];
while isempty(run_first_distances) || any(diff(run_first_distances)==0)
    run_first_distances = [randperm(3) randperm(3)]; % first distances in each of the 6 runs
end

% RULE: Block distances never repeat. Each combination of A-follows-B
% occurs 1-2 times per run. With 12 blocks, there will always be 1
% combination left over.
%
% Min number of each combination allowed in each run, distance x distance:
%   0     1     1
%   1     0     1
%   1     1     0
%
% Max number of each combination allowed in each run, distance x distance:
%   0     2     2
%   2     0     2
%   2     2     0
%
min_distance_combo_allowed = ones(3,3) * 1; % limit for each run
min_distance_combo_allowed(eye(3)==1) = 0;
max_distance_combo_allowed = ones(3,3) * 2; % limit for each run
max_distance_combo_allowed(eye(3)==1) = 0;

% RULE: The 1 left-over mentioned above must balanace across the 6 runs
%
% Number of each imbalance allowed across runs, distance x distance:
%   0     1     1
%   1     0     1
%   1     1     0
%
distance_combo_imbalance_remaining = ones(3,3) * 1;
distance_combo_imbalance_remaining(eye(3)==1) = 0;

% coordinates for rapid order calculation
[xs,ys] = meshgrid(1:3); % 3 distances x 3 distances

% initialize
run_distances_block = nan(12,6); % 12 blocks x 6 runs

% generate runs using a simple random search method
for run = 1:6
    % loop until success (cannot get stuck, there is always a common solution to this search)
    while 1
        % RULE: each distance occurs exactly twice in each half
        distances = [mod(randperm(6),3)+1 mod(randperm(6),3)+1];

        % must begin with pre-determined distance, else try again
        if distances(1) ~= run_first_distances(run)
            continue
        end

        % must no contain repeats, else try again
        if any(diff(distances)==0)
            continue
        end
    
        % calculate order matrix
        order_matrix = arrayfun(@(x,y) sum((distances(1:(end-1))==x)&(distances(2:end)==y)), xs, ys);
    
        % if any combinations less often than min, try again
        if any(order_matrix(:) < min_distance_combo_allowed(:))
            continue
        end

        % if any combinations more often than max, try again
        if any(order_matrix(:) > max_distance_combo_allowed(:))
            continue
        end

        % if the cell with a 1 has already had the imbalance, try again
        imbalance = (order_matrix==1);
        if ~distance_combo_imbalance_remaining(imbalance)
            continue
        end

        % success! note the order and cell with the imbalance
        run_distances_block(:,run) = distances;
        distance_combo_imbalance_remaining(imbalance) = distance_combo_imbalance_remaining(imbalance) - 1;
        break
    end
end

% check: there should not be any distance_combo_imbalance_remaining
if any(distance_combo_imbalance_remaining(:))
    error
end

% check: should not be any NaN in run_distances_block
if any(isnan(run_distances_block(:)))
    error
end


%% Part 2: Order of face sizes
% This process will always find a solution without needing to restart

fprintf("\tGenerating face size trial orders...\n");

% RULE: Independently for each distance, each size follows itself exactly 2
% times per run and follows each other size 2-3 times per run. This
% imbalance is controlled across runs by the next rule.
%
% Min number of each combination allowed in each run, size x size:
%   2     2     2
%   2     2     2
%   2     2     2
%
% Max number of each combination allowed in each run, size x size:
%   2     3     3
%   3     2     3
%   3     3     2
%
min_size_combo_allowed = ones(3,3) * 2; % min for each run
max_size_combo_allowed = ones(3,3) * 3; % max for each run
max_size_combo_allowed(eye(3)==1) = 2;

% RULE: In each run, 2 combinations will occur 3 times instead of 2. This
% imbalance will occur in each off-diagonal cell exactly twice across the 
% runs (and never in the diagonal cells).
%
% Number of each imbalance allowed across runs, size x size:
%   0     2     2
%   2     0     2
%   2     2     0
%
size_combo_imbalance_remaining = ones(3,3) * 2;
size_combo_imbalance_remaining(eye(3)==1) = 0;

% coordinates for rapid order calculation
[xs,ys] = meshgrid(1:3); % 3 sizes x 3 sizes

% initialize
run_sizes_trial = nan(6,12,6); % 6 trials x 12 blocks x 6 runs

% generate...
for distance = 1:3
    % initialize imbalance_remaining for this distance
    imbalance_remaining = size_combo_imbalance_remaining;

    % generate runs using a simple random search method
    for run = 1:6
        % find the blocks for this distance in this run
        ind_blocks = find(run_distances_block(:,run) == distance);

        % loop until success (cannot get stuck, there is always a semi-common solution to this search)
        % note that it could get stuck if imbalance_remaining was instead satisfied in a random order
        while 1
            % RULE: each size occurs twice per block (3 sizes x 2 reps =
            % 6 trials per block), recall: there are 4 blocks per distance
            sizes = [mod(randperm(6),3)+1 nan mod(randperm(6),3)+1 nan mod(randperm(6),3)+1 nan mod(randperm(6),3)+1]; % NaN are inserted between block as these trials do not follow each other
            
            % calculate order matrix
            order_matrix = arrayfun(@(x,y) sum((sizes(1:(end-1))==x)&(sizes(2:end)==y)), xs, ys);

            % if any combinations less often than min, try again
            if any(order_matrix(:) < min_size_combo_allowed(:))
                continue
            end
    
            % if any combinations more often than max, try again
            if any(order_matrix(:) > max_size_combo_allowed(:))
                continue
            end

            % if the cells with a 3 have already had the imbalance twice, try again
            imbalance = (order_matrix==3);
            if any(imbalance_remaining(imbalance) < max(imbalance_remaining(:))) % must solve the ones with 2 remaining before moving on to the 1s to avoid getting stuck and looping forever
                continue
            end

            % success! note the order and cell with the imbalance
            imbalance_remaining(imbalance) = imbalance_remaining(imbalance) - 1;
            run_sizes_trial(:, ind_blocks(1), run) = sizes(1:6);   % block 1
            run_sizes_trial(:, ind_blocks(2), run) = sizes(8:13);  % block 2
            run_sizes_trial(:, ind_blocks(3), run) = sizes(15:20); % block 3
            run_sizes_trial(:, ind_blocks(4), run) = sizes(22:27); % block 4
            break
        end
    end

    % check: there should not be any imbalance_remaining
    if any(imbalance_remaining(:))
        error
    end
end

% check: should not be any NaN in run_sizes_trial
if any(isnan(run_sizes_trial(:)))
    error
end

%% Part 3: Combine 3 sizes and 3 distances into 9 condition IDs

fprintf("\tGenerating lookup for 9 condition IDs (3 sizes x 3 distances)...\n");

% match 3 sizes x 3 distance to IDs 1-9
%
% size (rows) x distance (cols):
%   1     4     7
%   2     5     8
%   3     6     9
%
cond_IDs = nan(3,3);
cond_IDs(:) = 1:9;

% convert [12 block x 6 run] block distances into [6 trial x 12 block x 6 run] trial distances for per-trial lookup
run_distances_trial = repmat(reshape(run_distances_block,[1 12 6]), [6 1 1]);

% lookup condition IDs
run_IDs_trial = arrayfun(@(size,distance) cond_IDs(size,distance), run_sizes_trial, run_distances_trial);


%% Part 4: Order of face IDs

fprintf("\tGenerating face trial orders (longest step)...\n");

% RULE: In each run, each of 18 faces will occur during each of 9 condition
% IDs (3 size x 3 distance) 0-1 times.
%
min_face_per_ID_run = 0;
max_face_per_ID_run = 1;

% RULE: Across all runs, each of 18 faces will occur during each of 9
% condition IDs (3 size x 3 distance) 1-3 times.
%
min_face_per_ID_total = 1; % 2 would be better but solution is rare or maybe impossible when combined with other rules
max_face_per_ID_total = 3;

% RULE: In each run, each of 18 faces will (a) not follow itself and (b)
% follow each other faces 0-1 times.
%
% Min number of each combination allowed in each run, face x face:
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
%
% Max number of each combination allowed in each run, face x face:
% 0     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
% 1     0     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
% 1     1     0     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
% 1     1     1     0     1     1     1     1     1     1     1     1     1     1     1     1     1     1
% 1     1     1     1     0     1     1     1     1     1     1     1     1     1     1     1     1     1
% 1     1     1     1     1     0     1     1     1     1     1     1     1     1     1     1     1     1
% 1     1     1     1     1     1     0     1     1     1     1     1     1     1     1     1     1     1
% 1     1     1     1     1     1     1     0     1     1     1     1     1     1     1     1     1     1
% 1     1     1     1     1     1     1     1     0     1     1     1     1     1     1     1     1     1
% 1     1     1     1     1     1     1     1     1     0     1     1     1     1     1     1     1     1
% 1     1     1     1     1     1     1     1     1     1     0     1     1     1     1     1     1     1
% 1     1     1     1     1     1     1     1     1     1     1     0     1     1     1     1     1     1
% 1     1     1     1     1     1     1     1     1     1     1     1     0     1     1     1     1     1
% 1     1     1     1     1     1     1     1     1     1     1     1     1     0     1     1     1     1
% 1     1     1     1     1     1     1     1     1     1     1     1     1     1     0     1     1     1
% 1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     0     1     1
% 1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     0     1
% 1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     0
%
min_face_order_run = zeros(18,18);
max_face_order_run = ones(18,18);
max_face_order_run(eye(18)==1) = 0;

% RULE: Across all runs, each of 18 faces will (a) not follow itself and 
% (b) follow each other faces 1-2 times.
%
% Min number of each combination allowed total, face x face:
% 0     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
% 1     0     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
% 1     1     0     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
% 1     1     1     0     1     1     1     1     1     1     1     1     1     1     1     1     1     1
% 1     1     1     1     0     1     1     1     1     1     1     1     1     1     1     1     1     1
% 1     1     1     1     1     0     1     1     1     1     1     1     1     1     1     1     1     1
% 1     1     1     1     1     1     0     1     1     1     1     1     1     1     1     1     1     1
% 1     1     1     1     1     1     1     0     1     1     1     1     1     1     1     1     1     1
% 1     1     1     1     1     1     1     1     0     1     1     1     1     1     1     1     1     1
% 1     1     1     1     1     1     1     1     1     0     1     1     1     1     1     1     1     1
% 1     1     1     1     1     1     1     1     1     1     0     1     1     1     1     1     1     1
% 1     1     1     1     1     1     1     1     1     1     1     0     1     1     1     1     1     1
% 1     1     1     1     1     1     1     1     1     1     1     1     0     1     1     1     1     1
% 1     1     1     1     1     1     1     1     1     1     1     1     1     0     1     1     1     1
% 1     1     1     1     1     1     1     1     1     1     1     1     1     1     0     1     1     1
% 1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     0     1     1
% 1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     0     1
% 1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     0
%
% Max number of each combination allowed total, face x face:
% 0     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2
% 2     0     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2
% 2     2     0     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2
% 2     2     2     0     2     2     2     2     2     2     2     2     2     2     2     2     2     2
% 2     2     2     2     0     2     2     2     2     2     2     2     2     2     2     2     2     2
% 2     2     2     2     2     0     2     2     2     2     2     2     2     2     2     2     2     2
% 2     2     2     2     2     2     0     2     2     2     2     2     2     2     2     2     2     2
% 2     2     2     2     2     2     2     0     2     2     2     2     2     2     2     2     2     2
% 2     2     2     2     2     2     2     2     0     2     2     2     2     2     2     2     2     2
% 2     2     2     2     2     2     2     2     2     0     2     2     2     2     2     2     2     2
% 2     2     2     2     2     2     2     2     2     2     0     2     2     2     2     2     2     2
% 2     2     2     2     2     2     2     2     2     2     2     0     2     2     2     2     2     2
% 2     2     2     2     2     2     2     2     2     2     2     2     0     2     2     2     2     2
% 2     2     2     2     2     2     2     2     2     2     2     2     2     0     2     2     2     2
% 2     2     2     2     2     2     2     2     2     2     2     2     2     2     0     2     2     2
% 2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     0     2     2
% 2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     0     2
% 2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     0
%
min_face_order_total = ones(18,18) * 1;
min_face_order_total(eye(18)==1) = 0;
max_face_order_total = ones(18,18) * 2;
max_face_order_total(eye(18)==1) = 0;

% generate with heuristic search...
%   i.e., randomly selections from the best options at each step and resets whenever it gets stuck
%   WARNING: editing this section may result in a loop that never ends, which will likely freeze up the computer
%   NOTE: solutions to these rules are quite rare and may require multiple minutes of searching
max_attempts_per_run = 100; % after this many attempts on a single run, the entire process will reset
while 1

% initialize
face_per_ID_total = zeros(18,9); % 18 face x 9 condition ID (must be between min_face_per_ID_total and max_face_per_ID_total)
face_order_total = zeros(18,18); % 18 face x 18 face        (must be between min_face_order_total  and max_face_order_total)
valid_all = true;
run_faces_trial = nan(6,12,6); % 6 trials x 12 blocks x 6 runs

for run = 1:6
    c=0;
    while 1
        c=c+1;
        if c>max_attempts_per_run
            % probably stuck in a bad order, reset
            valid_all = false;
            break;
        end

        % init
        face_per_ID_run = zeros(18,9); % 18 face x 9 condition ID (must be between min_face_per_ID_run and max_face_per_ID_run)
        face_order_run = zeros(18,18); % 18 face x 18 face        (must be between min_face_order_run  and max_face_order_run)
        valid_run = true;
        face_per_ID_total_temp = face_per_ID_total;
        face_order_total_temp = face_order_total;

        % default to no faces available, will trigger filling with 18 faces
        faces_available = false(18,1);
    
        for block = 1:12
            % RULE: each face will occur exactly once in each quarter of a 
            % run. This this is implemented by refilling the 18 faces once
            % they are all used
            if ~any(faces_available)
                faces_available = true(18,1);
            end

            face_prior = nan;
            for trial = 1:6
                % which cond ID is this trial?
                ID = run_IDs_trial(trial, block, run);

                % initialize to available faces
                face_valid = faces_available;

                % run level: which faces can be paired with this ID?
                face_valid = face_valid & (face_per_ID_run(:,ID) < max_face_per_ID_run);

                % total: which faces can be faired with this ID?
                face_valid = face_valid & (face_per_ID_total_temp(:,ID) < max_face_per_ID_total);

                % run level: which faces can follow the prior face?
                if ~isnan(face_prior)
                    face_valid = face_valid & (face_order_run(:,face_prior) < max_face_order_run(:,face_prior));
                end

                % total: which faces can follow the prior face?
                if ~isnan(face_prior)
                    % this bit of code prioratizes options that fulfill face_order_total_temp==0
                    x = face_order_total_temp(:,face_prior);
                    x(face_prior) = inf;
                    x(~face_valid) = inf;
                    if any(x==0)
                        face_valid = face_valid & (face_order_total_temp(:,face_prior) == 0); 
                    end

                    % applies the standard rules for maximums
                    face_valid = face_valid & (face_order_total_temp(:,face_prior) < max_face_order_total(:,face_prior));
                end

                % if there aren't any options left, restart the run
                if ~any(face_valid)
                    valid_run = false;
                    break; % if invalid, stop checking trials
                end

                % randomly select one of the options
                ind_options = find(face_valid);
                ind_select = randperm(length(ind_options),1);
                face_select = ind_options(ind_select);

                % apply selection
                faces_available(face_select) = false;
                run_faces_trial(trial, block, run) = face_select;
                face_per_ID_run(face_select, ID) = face_per_ID_run(face_select, ID) + 1;
                face_per_ID_total_temp(face_select, ID) = face_per_ID_total_temp(face_select, ID) + 1;
                if ~isnan(face_prior)
                    face_order_run(face_select,face_prior) = face_order_run(face_select,face_prior) + 1;
                    face_order_total_temp(face_select,face_prior) = face_order_total_temp(face_select,face_prior) + 1;
                end
                face_prior = face_select;
            end

            % if invalid, stop checking runs
            if ~valid_run
                break;
            end
        end
        % if invalid, loop back to try the run again
        if ~valid_run
            continue
        end

        % check min_face_per_ID_run
        if any(face_per_ID_run(:) < min_face_per_ID_run)
            continue
        end

        % check min_face_order_run
        if any(face_order_run < min_face_order_run)
            continue
        end
    
        % successful run
        break

    end

    % if overall invalid, stop checking runs
    if ~valid_all
        break
    end

    % update running totals
    face_per_ID_total = face_per_ID_total_temp;
    face_order_total = face_order_total_temp;
end

% if overall invalid, loop back and try again
if ~valid_all
    continue
end

% check min_face_per_ID_total
if any(face_per_ID_total(:) < min_face_per_ID_total)
    continue
end

% check min_face_order_total
if any(face_order_total < min_face_order_total)
    continue
end

% successful set of runs
break
end

%% Part 5: Combine everything into an order spreadsheet

fprintf("\tGenerating order table...\n");

% Reminders:
% block distances:  run_distances_block     12 block x 6 run
% trial sizes:      run_sizes_trial         6 trial x 12 block x 6 run
% trial faces:      run_faces_trial         6 trial x 12 block x 6 run
% trial IDs:        run_IDs_trial           6 trial x 12 block x 6 run

% for each run...
for run = 1:6

% initialize table
tbl = table(Size=[146 11], ...
            VariableNames = ["Trial" "Condition" "Duration_Seconds" "Filename_Left" "Filename_Right" "Is_Repeat" "Block" "Distance" "Face" "Size" "ConditionID"], ...
            VariableTypes = ["double" "string" "double" "string" "string" "logical" "double" "string" "string" "string" "double"]);
row = 0;

% add each event...

% initial baseline (-1)
row = row + 1;
tbl.Trial(row) = 0;
tbl.Condition(row) = "NULL_Baseline";
tbl.Duration_Seconds(row) = durations.baseline_initial;
tbl.Filename_Left(row) = stims.image_prefix_baseline + "_" + stims.baseline_distance + "_" + stims.image_suffix_left_eye + stims.image_filetype;
tbl.Filename_Right(row) = stims.image_prefix_baseline + "_" + stims.baseline_distance + "_" + stims.image_suffix_right_eye + stims.image_filetype;
tbl.Is_Repeat(row) = false;
tbl.Block(row) = 0;
tbl.Distance(row) = stims.baseline_distance;
tbl.Face(row) = "NULL";
tbl.Size(row) = "NULL";
tbl.ConditionID(row) = "-1";

% blocks...
trial_counter_total = 0;
for block = 1:12
    distance = stims.distances(run_distances_block(block, run));

    % vergence change (-2)...
    row = row + 1;
    tbl.Trial(row) = 0;
    tbl.Condition(row) = "NULL_Change";
    tbl.Duration_Seconds(row) = durations.baseline_vergence_change;
    tbl.Filename_Left(row) = stims.image_prefix_vergence_change + "_" + distance + "_" + stims.image_suffix_left_eye + stims.image_filetype;
    tbl.Filename_Right(row) = stims.image_prefix_vergence_change + "_" + distance + "_" + stims.image_suffix_right_eye + stims.image_filetype;
    tbl.Is_Repeat(row) = false;
    tbl.Block(row) = block;
    tbl.Distance(row) = distance;
    tbl.Face(row) = "NULL";
    tbl.Size(row) = "NULL";
    tbl.ConditionID(row) = "-2";

    % trials...
    for trial = 1:6
        size = stims.sizes(run_sizes_trial(trial,block,run));
        face = stims.faces(run_faces_trial(trial,block,run));
        condID = run_IDs_trial(trial,block,run);
        trial_counter_total = trial_counter_total + 1;

        % stim...
        row = row + 1;
        tbl.Trial(row) = trial_counter_total;
        tbl.Condition(row) = face + "_" + distance + "_" + size;
        tbl.Duration_Seconds(row) = durations.stim;
        tbl.Filename_Left(row) = tbl.Condition(row) + "_" + stims.image_suffix_left_eye + stims.image_filetype;
        tbl.Filename_Right(row) = tbl.Condition(row) + "_" + stims.image_suffix_right_eye + stims.image_filetype;
        tbl.Is_Repeat(row) = false;
        tbl.Block(row) = block;
        tbl.Distance(row) = distance;
        tbl.Face(row) = face;
        tbl.Size(row) = size;
        tbl.ConditionID(row) = condID;

        % ITI (-3)...
        if trial < 6
            row = row + 1;
            tbl.Trial(row) = 0;
            tbl.Condition(row) = "NULL_ITI";
            tbl.Duration_Seconds(row) = durations.ITI;
            tbl.Filename_Left(row) = stims.image_prefix_ITI + "_" + distance + "_" + stims.image_suffix_left_eye + stims.image_filetype;
            tbl.Filename_Right(row) = stims.image_prefix_ITI + "_" + distance + "_" + stims.image_suffix_right_eye + stims.image_filetype;
            tbl.Is_Repeat(row) = false;
            tbl.Block(row) = block;
            tbl.Distance(row) = distance;
            tbl.Face(row) = "NULL";
            tbl.Size(row) = "NULL";
            tbl.ConditionID(row) = -3;
        end
    end
end

% final baseline (-4)
row = row + 1;
tbl.Trial(row) = 0;
tbl.Condition(row) = "NULL_Baseline";
tbl.Duration_Seconds(row) = durations.baseline_final;
tbl.Filename_Left(row) = stims.image_prefix_baseline + "_" + stims.baseline_distance + "_" + stims.image_suffix_left_eye + stims.image_filetype;
tbl.Filename_Right(row) = stims.image_prefix_baseline + "_" + stims.baseline_distance + "_" + stims.image_suffix_right_eye + stims.image_filetype;
tbl.Is_Repeat(row) = false;
tbl.Block(row) = 0;
tbl.Distance(row) = stims.baseline_distance;
tbl.Face(row) = "NULL";
tbl.Size(row) = "NULL";
tbl.ConditionID(row) = "-4";

% expected size?
if row ~= height(tbl)
    error
end

% save
fp = fol + sprintf("PAR%02d_RUN%02d.xlsx", par, run);
fprintf("\tSaving: %s\n", fp);
if exist(fp, "file")
    delete(fp)
end
writetable(tbl, fp);

end


%% end of participant loop
end
%% done
disp Done!