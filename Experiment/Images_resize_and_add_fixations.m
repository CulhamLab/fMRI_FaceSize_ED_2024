fol = ".\Images\";
fol_fixation = fol + "Fixation_To_Add\";
fol_without_fixation = fol + "Stims_Without_Fixation\";

threshold = 20;

img_width = 1920;
img_height = 1080;

for eye = ["l" "r"]
    for distance = [51 64 80]
        name = sprintf("Fixation_D%02d_%s.png", distance, eye);
        fix = imread(fol_fixation + name);
        fgd = repmat(any(fix <= threshold,3),[1 1 3]);

        % resize
        fix_rs = imresize(fix, [img_height img_width]);
        
        % save
        imwrite(fix_rs, fol + name);

        list = dir(fol_without_fixation + sprintf("*_D%02d_*%s.png", distance, eye));
        for file = list(:)'
            % add fixation
            img = imread([file.folder filesep file.name]);
            img(fgd) = fix(fgd);

            % resize
            img = imresize(img, [img_height img_width]);
            
            % save
            imwrite(img, fol + file.name);
        end
    end
end

