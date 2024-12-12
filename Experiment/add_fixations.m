fol = ".\Images\";
fol_without_fixation = fol + "WithoutFixation\";

for eye = ["L" "R"]
    for distance = [51 64 80]
        fix = imread(fol + sprintf("Fixation_D%02d_%s.png", distance, eye));
        fgd = repmat(any(fix,3),[1 1 3]);

        list = dir(fol_without_fixation + sprintf("Face*_D%02d_S*_%s.png", distance, eye));
        for file = list(:)'
            img = imread([file.folder filesep file.name]);
            img(fgd) = fix(fgd);
            imwrite(img, fol + file.name);
        end
    end
end

