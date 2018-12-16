function colorGradient = getAlphaGradient(color1, color2, n)

    % initialize
    colorGradient = nan(n, 3);
    
    % compute gradient
    colorGradient(:,1) = linspace(color1(1),color2(1),n);
    colorGradient(:,2) = linspace(color1(2),color2(2),n);
    colorGradient(:,3) = linspace(color1(3),color2(3),n);
    
%     colorGradient = colorGradient(fliplr(1:size(colorGradient,1)),:);
end