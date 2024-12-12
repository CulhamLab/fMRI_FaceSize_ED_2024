%%
fol = ".\Images\";
if ~exist(fol, "dir"), mkdir(fol); end
fol_fixation = fol + "Fixation_To_Add\";
if ~exist(fol_fixation, "dir"), mkdir(fol_fixation); end
fol_without_fixation = fol + "Stims_Without_Fixation\";
if ~exist(fol_without_fixation, "dir"), mkdir(fol_without_fixation); end

%%
img_width = 1920;
img_height = 1080;

font_size = 30;
font_colour = [1 1 1];
colour_left = [0 0 1];
colour_right = [1 0 0];

eyes = ["l" "r"];
distances = [51 64 80];
sizes = [20 25 31];

for eye = ["l" "r"]
    for distance = distances
        % Fixation
        y = (0.1 + ((find(distances==distance)-1)/2*0.8)) * img_height;
        img = ones(img_height,img_width,3,"uint8") * 255;
        switch eye
            case "l"
                img = insertText(img, [(0.4 * img_width) y], "Fixation " + distance, FontSize=font_size, TextBoxColor=[1 1 1], FontColor=[0 0 0]);
            case "r"
                img = insertText(img, [(0.5 * img_width) y], "Fixation " + distance, FontSize=font_size, TextBoxColor=[1 1 1], FontColor=[0 0 0]);
        end
        imwrite(img, fol_fixation + sprintf("Fixation_D%02d_%s.png", distance, eye));

        % Ball
        img = ones(img_height,img_width,3,"uint8") * 128;
        switch eye
            case "l"
                img = insertText(img, [(0.1 * img_width) y], "Ball " + distance, FontSize=fsz, TextBoxColor=colour_left, FontColor=font_colour);
            case "r"
                img = insertText(img, [(0.7 * img_width) y], "Ball " + distance, FontSize=fsz, TextBoxColor=colour_right, FontColor=font_colour);
        end
        imwrite(img, fol_without_fixation + sprintf("Ball_D%02d_%s.png", distance, eye));

        % Faces
        for sz = sizes
            name = sprintf("Face%02d_D%02d_S%02d_%s.png", 1, distance, sz, eye);
            img = ones(img_height,img_width,3,"uint8") * 128;
            fsz = ceil(font_size * 1.5 * (sz / max(sizes)));
            switch eye
                case "l"
                    img = insertText(img, [(0.1 * img_width) y], name, FontSize=fsz, TextBoxColor=colour_left, FontColor=font_colour);
                case "r"
                    img = insertText(img, [(0.7 * img_width) y], name, FontSize=fsz, TextBoxColor=colour_right, FontColor=font_colour);
            end
            imwrite(img, fol_without_fixation + name);

        end
    end
end