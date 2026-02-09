function [cmap,inColorLocations] = makeUniformColormap(rgb_in,N)
% Takes a colormap and adjusts the spacing between colors to make it
% roughly perceptually uniform in LAB space.
%
% (c) Blox Bloxham 2026

if any(rgb_in > 1,'all') && all(rgb_in >= 0,'all') && all(rgb_in <= 255,'all')
    warning('Assuming RGB input supplied on [0,255] scale.  Output will be [0,1]-scale RGB.')
    lab_in = rgb_in/255;
elseif any(rgb_in > 1,'all') || any(rgb_in < 0,'all')
    warning('Assuming input supplied as CIELAB values.  Output will be RGB.')
    lab_in = rgb_in;
else
    lab_in = rgb2lab(rgb_in);
end

nonuniform = rgb2lab(min(max(lab2rgb([...
    interp1(linspace(0,1,size(lab_in,1)),lab_in(:,1),linspace(0,1,N));
    interp1(linspace(0,1,size(lab_in,1)),lab_in(:,2),linspace(0,1,N));
    interp1(linspace(0,1,size(lab_in,1)),lab_in(:,3),linspace(0,1,N))]'),...
    0),1));

dist = [0;sqrt(sum((nonuniform(2:end,:) - nonuniform(1:end-1,:)).^2,2))];

cumdist = cumsum(dist);

cmap = min(max(lab2rgb([...
    interp1(cumdist,nonuniform(:,1),linspace(0,cumdist(end),N));
    interp1(cumdist,nonuniform(:,2),linspace(0,cumdist(end),N));
    interp1(cumdist,nonuniform(:,3),linspace(0,cumdist(end),N))]'),0),1);

if nargout > 1
    inColorLocations = zeros(size(rgb_in,1),1);
    for i = 1:size(rgb_in,1)
        [~,inColorLocations(i)] = min(sum((cmap - rgb_in(i,:)).^2,2));
    end
end

end

