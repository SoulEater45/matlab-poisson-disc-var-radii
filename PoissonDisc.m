clear;
close all;

dimensions = [640, 360];
dimensions = [256, 256, 256];
n = length(dimensions);
r = 30;
radii = [10 20];
k = 50;
%w = r / sqrt(n);
w = 2*min(radii) / sqrt(n);
o = 5;

vis = true;

if all([vis any(n==[2,3])])
    figure;
    set(gcf, 'WindowState', 'maximized');
    cmap = jet(length(radii));
    %pause(5);
    % if vis && n == 2
    % 
    % elseif vis && n==3 && vis
    % 
    % end
end
%% Step 0
bgGrid = cell(ceil(dimensions / w));
bgGridRadii = zeros(ceil(dimensions / w));

%% Step 1
pos = rand(1,n) .* dimensions;
p = ceil(pos/w);
p = eval(sprintf('sub2ind(size(bgGrid), %s);', join("p(" + string(1:length(p)) + ")", ",")));
bgGrid{p} = pos;
bgGridRadii(p) = radii(randi(length(radii)));
active = cell(1,1);
active{1} = pos;
activeRadii = bgGridRadii(p);

if all([vis any(n==[2,3])])
    axis equal;
    if n == 2
        viscircles(bgGrid{p}, bgGridRadii(p));
        rectangle('Position', [bgGrid{p}, 3*bgGridRadii(p), 3*bgGridRadii(p)]-bgGridRadii(p), ...
            'Curvature', [1 1], ...
            'FaceColor', cmap(bgGridRadii(p) == radii, :));
        axis([0, dimensions(1), 0, dimensions(2)]);
        %xticks(0:w:dimensions(1));
        %yticks(0:w:dimensions(2));
        %set(gca,'xticklabel',[]);
        %set(gca,'yticklabel',[]);
        %grid on;
    elseif n == 3
        [xS, yS, zS] = sphere;
        surf(...
            bgGrid{p}(1)+xS.*bgGridRadii(p), ...  % Shift and scale x data
            bgGrid{p}(2)+yS.*bgGridRadii(p), ...  % Shift and scale y data
            bgGrid{p}(2)+zS.*bgGridRadii(p), ...  % Shift and scale z data
            'EdgeColor', 'none', ...
            'FaceColor', cmap(bgGridRadii(p) == radii, :));
        axis([0, dimensions(1), 0, dimensions(2), 0, dimensions(3)]);
        %xticks(0:w:dimensions(1));
        %yticks(0:w:dimensions(2));
        zticks(0:w:dimensions(3));
        %set(gca,'xticklabel',[]);
        %set(gca,'yticklabel',[]);
        set(gca,'zticklabel',[]);
        %grid on;
        %hold on;
    end
    xticks(0:w:dimensions(1));
    yticks(0:w:dimensions(2));
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    grid on;
    hold on;
end

%% Step 2
while ~isempty(active)
    activeIdx = randi(length(active));
    activePoint = active{activeIdx};
    activeRadius = activeRadii(activeIdx);
    
    % TODO: generate k samples and do the maths parallel
    found = false;
    for j = 1:k
        % Sample is generated close to the active point with a distance
        % of r*(1+rand) which is between r and 2r
        sampleRadius = radii(randi(length(radii)));
        sample = activePoint + HyperSphere(2*pi*rand(n-1,1), sampleRadius + activeRadius + min(sampleRadius,activeRadius)*rand)';
        
        if any([sample < 0, sample > dimensions])
            continue;
        end
        
        % if n == 2
        %     hold on;
        %     viscircles(sample, r, 'Color', 'g');
        %     plot(sample(1), sample(2), 'g.', 'MarkerSize', 10);
        %     hold off;
        % end
        
        d = ceil((max(radii)+sampleRadius+o)/w);
        p = ceil(sample / w);
        p = p - d + (linspace(0,2*d,2*d+1) .* ones(n, 2*d+1))';
        mgrid = zeros([n (2*d+1) * ones(1,n)]);
        eval(sprintf('[%s] = ndgrid(%s);', join("mgrid(" + string(1:n) + repmat(',:', 1, n) + ")", ","), join("p(:," + string(1:n) + ")", ",")));
        mgrid = reshape(mgrid, n, []);
        %mgrid = mgrid(:,~any([mgrid<1;mgrid>size(bgGrid)';all(mgrid == ceil(sample/w)',1)],1));
        mgrid = mgrid(:,~any([mgrid<1;mgrid>size(bgGrid)'],1));
        eval(sprintf('p = sub2ind(size(bgGrid),%s);', join(" mgrid(" + string(1:n) + repmat(',:', 1, n) + ")",",")));
        
        collision = false;
        if ~isempty(cell2mat(bgGrid(p)'))
            collision = any(sum((cell2mat(bgGrid(p)')-sample).^2, 2) < (rand*o + sampleRadius + bgGridRadii(p)).^2);
        end
        
        if ~collision
            p = ceil(sample/w);
            p = eval(sprintf('sub2ind(size(bgGrid), %s);', join("p(" + string(1:length(p)) + ")", ",")));
            bgGrid{p} = sample;
            active{end+1} = sample;
            bgGridRadii(p) = sampleRadius;
            activeRadii(end+1) = sampleRadius;
            found = true;
            break;
        end
        
        if vis && n == 3
            v = get(gca, 'View');
            v(1) = v(1) + 1/k;
            view(v);
        end
    end
    
    if ~found
        active(activeIdx) = [];
        activeRadii(activeIdx) = [];
    end
    
    if all([found vis any(n==[2,3])])
        if n==2
            %viscircles(cell2mat(active'), ones(length(active), 1)*r/2);
            %viscircles(sample, sampleRadius);
            rectangle('Position', [sample, 3*sampleRadius, 3*sampleRadius]-sampleRadius, ...
                'Curvature', [1 1], ...
                'FaceColor', cmap(sampleRadius == radii, :));
            drawnow;
        elseif n==3
            surf(...
                sample(1)+xS.*sampleRadius, ...  % Shift and scale x data
                sample(2)+yS.*sampleRadius, ...  % Shift and scale y data
                sample(3)+zS.*sampleRadius, ...  % Shift and scale z data
                'EdgeColor', 'none', ...
                'FaceColor', cmap(sampleRadius == radii, :));
            axis equal;
            axis([0, dimensions(1), 0, dimensions(2), 0, dimensions(3)]);
        end
        pause(0.01);
    end
end