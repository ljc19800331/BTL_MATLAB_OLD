%Input: An RGB image
%Output: [x,y] coordinates of points corrresponding to the selected color;
%        image_dilated - a binary image showing the points as 1's

function [xy_coordinates, image_dilated] = process_tattoo_image(filename)
%% Process tattoo image
    myTattoo = imread(filename);
    myTattooColor = myTattoo;
    h1 = figure; imshow(myTattooColor); title('Original in Grayscale.');

    disp('Select a color for search and press enter.')

    [x,y] = getpts();

    target_color = myTattooColor(uint16(y),uint16(x),:);

    target_color = reshape(target_color, [1,1,3]);

    difference = single(myTattooColor) - single(target_color);

    error = sqrt(sum(difference.^2,3));

    filterImage = error < 10;

    h0 = figure;
    hold on;
    slider_h = uicontrol('style','slider');

%% Clean it up some.
    while(ishandle(h0))
        imshow(~filterImage); title('Determine proper threshold...');
        try
            value = get(slider_h,'value')*150;
        catch
            disp('done');
        end
        %filterImage = sqrt(sum((myTattooColor - color_val).^2,3)) <= value;
        filtered = error <= value;
        strel = [0 1 0; 1 1 1; 0 1 0];
        strel = [1];
        filterImage = imerode(filtered,strel);
        disp(value)
        pause(0.01)
    end

%% clean up the image by hand.
    disp('Clean image.');
    h1 = figure; title('Clean the image...')
    while(ishandle(h1))
        imshow(filterImage);
        try
            disp('Select region to delete...');
            rect = getrect(h1);
        catch
            disp('caught');
        end

        rect = floor(rect);
        xmin = max(rect(2),1);
        xmax = min(rect(2) + rect(4),size(filterImage,1));
        ymin = max(rect(1),1);
        ymax = min(rect(1) + rect(3),size(filterImage,2));
        filterImage(xmin:xmax,ymin:ymax) = 0;
    end

    filter_new = filterImage;
    h2 = figure; imshow(filter_new); title('Filter New');

%% Get image coordinates
    [M1, M2] = ind2sub(size(filter_new), find(filter_new == 1));
    hold on; scatter(M2, M1, 1);
    xy_coordinates = [M1, M2];

    image_dilated = filter_new;
end
