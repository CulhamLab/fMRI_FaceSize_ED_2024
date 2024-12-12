%%
fol = ".\Images\";
if ~exist(fol, "dir"), mkdir(fol); end
fol_without_fixation = fol + "WithoutFixation\";
if ~exist(fol_without_fixation, "dir"), mkdir(fol_without_fixation); end

%%
img_width = 1920;
img_height = 1080;

font_size = 30;
font_colour = [1 1 1];
colour_left = [0 0 1];
colour_right = [1 0 0];

eyes = ["L" "R"];
distances = [51 64 80];
sizes = [20 25 31];

for eye = ["L" "R"]
    for distance = distances
        % Fixation
        y = (0.1 + ((find(distances==distance)-1)/2*0.8)) * img_height;
        img = zeros(img_height,img_width,3,"uint8");
        switch eye
            case "L"
                img = insertText(img, [(0.4 * img_width) y], "Fixation " + distance, FontSize=font_size, TextBoxColor=colour_left, FontColor=font_colour);
            case "R"
                img = insertText(img, [(0.5 * img_width) y], "Fixation " + distance, FontSize=font_size, TextBoxColor=colour_right, FontColor=font_colour);
        end
        imwrite(img, fol + sprintf("Fixation_D%02d_%s.png", distance, eye));

        % Faces
        for sz = sizes
            name = sprintf("Face%02d_D%02d_S%02d_%s.png", 1, distance, sz, eye);
            img = zeros(img_height,img_width,3,"uint8");
            fsz = ceil(font_size * 1.5 * (sz / max(sizes)));
            switch eye
                case "L"
                    img = insertText(img, [(0.1 * img_width) y], name, FontSize=fsz, TextBoxColor=colour_left, FontColor=font_colour);
                case "R"
                    img = insertText(img, [(0.7 * img_width) y], name, FontSize=fsz, TextBoxColor=colour_right, FontColor=font_colour);
            end
            imwrite(img, fol_without_fixation + name);

        end
    end
end