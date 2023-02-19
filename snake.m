%%%%%%CONTROL VARIABLES%%%%%%%%%%
ts = 0.03; %time step
g = 1; %damping coeffiecent

%internal forces
a = 0.2; % weight for alpha - resistance to streching deformations
b = 0.4; % weight for beta - resistance to bending deformations

%external forces
T = 0.001; % treshold for image intensity
qConst = 50; % amplitude of inflation force
qloweringStep = 2;
stuckIterations = 4;
stuckDistance = 0.1;

p = 51; % weight of external image force (should be slightly larger than q
c = 1; % scale of potential

%slices vatiables
numOfIter = 20; % iteration limit
numOfCircNodes = 10; % how many nodes of the active contour are to be generated
circRadius = 2; % size of the starting circle described by the generated nodes
everyNchoose = 257; % after how many slices do you want to manually set the center again

%enable debug view
debug = false;
sliceCounter = 0;

%%%%%%%%%%%%SNAKE ALGORITHM%%%%%%%%%%%%%
%FILE
%clear the contents of file
fid = fopen('output.txt', 'w+');
fprintf(fid, '');
fclose(fid);
%open file for writing
fid = fopen('output.txt', 'a');

%read nrrd slices
A = nhdr_nrrd_read('fib1321.seg.nrrd', true);

[d1, d2, d3] = size(A.data);
I = A.data(:, :, 1);
I = im2double(I);
imagesc(I);
xc = 0; yc = 0;
jSum = 0;
for k=1:d3
    fprintf('processing slice %i\n', sliceCounter)
    sliceCounter = sliceCounter + 1;
    q = qConst; 
    j = 0;
    prevInterF = 1;
    I = A.data(:, :, k);
    I = im2double(I);
    
    if(mod(k-1,everyNchoose) == 0)
        % Show the image and select a point with the mouse
        h = gca;
        h.Visible = 'On';
        fprintf('right click on the centre of the object you want to segment\n')
        [x,y] = getpts;
        fprintf('center selected\n')
        xc = x(1);
        yc =y(1);
        tic;
    end  
    
    %automaticly chose points in a circle
    [x,y] = getCircPoints(xc, yc, numOfCircNodes, circRadius);
    N=[x(:) y(:)];

    %m is number of points
    m = size(N, 1);

    %potential of image
    P = gradient(imgaussfilt(I));
    Pnorm = P/norm(P);
    P = -c*Pnorm;
    %imshow(P);

    Nnew = zeros(m, 2);
    qM = ones(m, 1)*stuckIterations;  
    interFM = zeros(m, 1);
  
    while true
        j = j + 1;
        %go through all the nodes and update their position based on acting forces
        % calculate normals of nodes
        xt=N(:,1); yt=N(:,2);
        dy = gradient(yt);
        dx = gradient(xt);
        %normalize
        len=sqrt(dx.^2+dy.^2); %lahko tudi norm
        normals = [dy./len -dx./len];

        % reset figure
        imagesc(I);
        h = gca;
        h.Visible = 'On';

        for i=1:m
            %calculate indexes so the boundary condition holds true (x1 = xm)
            iPr = i-1; iPrPr = i-2;
            iNx = i+1; iNxNx = i+2;
            if(i == 1)
                iPr = m;
                iPrPr = m-1;
            elseif(i == 2)
                iPrPr = m;
            elseif(i == m)
                iNx = 1;
                iNxNx = 2;
            elseif(i == m-1)
                iNxNx = 1;
            end  
            xi = N(i,:);
            xiPr = N(iPr,:); xiPrPr = N(iPrPr,:);
            xiNx = N(iNx,:); xiNxNx = N(iNxNx,:);

            %internal forces
            alfa = 2*xi-xiPr-xiNx; %tensile forces
            alfaPrev = 2*xiPr-xiPrPr-xi;
            alfaNext = 2*xiNx-xi-xiNxNx;

            beta = 2*alfa - alfaPrev - alfaNext; %flexural force

            %external forces
            %interpolate image intensity and treshold it with binary function
            interI = interp2(I, xi(1), xi(2), 'linear');
            if interI >= T, interF = 1; else, interF = -1; end % tresholding
            interFM(i) = interF; % save interF to matrix for later node end condition checking
            rho = q*interF*normals(i,:); % inflation force pushes the model towards intensity edges

            %interpolate potential
            interP = interp2(P, xi(1), xi(2), 'linear');
            ft = p*(gradient(interP)); %image force to stop the contour at significant edges

            Nnew(i, :) = xi - (ts/g)*(a*alfa + b*beta - rho - ft);

            %preventing oscillation by lowering q if the node changed direction
            if(j ~= 1 && prevInterF < 0 && interF > 0)
                q = q - qloweringStep;
            elseif(j ~= 1 && prevInterF > 0 && interF < 0)
                q = q - qloweringStep;
            end
            prevInterF = interF;

            if debug == true
                hold on; plot(xi(1), xi(2), 'c*');
                text(xi(1), xi(2), num2str(i));
                quiver(xi(1), xi(2), normals(i,1)*10, normals(i,2)*10, 'c');
                hold on; plot(Nnew(i,1),Nnew(i,2), 'r*');
                text(Nnew(i,1),Nnew(i,2), strcat('n',num2str(i)));
            end
        end
        
        %early break conditions
        if(j ~= 1) %break if half of nodes are stuck
            mask1 = abs(N-Nnew) < stuckDistance;
            mask2 = mask1(:,1) + mask1(:,2);
            mask2 = mask2 == 2;

            qM = qM-mask2;
            numOfStuck = sum(qM < 0);
            if(numOfStuck >= m/2)
                break;
            end
        end

        if(q < 0) %break if energy is zero
           break;
        end

        if(j > numOfIter) %break if iteration limits for one slice is reached
            break;
        end

        N = Nnew;
    end
    
    jSum = jSum + j;
    
    %center of current 'polygon'
    xc = mean(N(:,1));
    yc = mean(N(:,2));
    if debug == true
        hold on; plot(xc,yc, 'b+');
    end

    % finish if all nodes are out of object
    if(sum(interFM < 0) == m)
        break
    end
    
    for i=1:m  
        fprintf(fid, '%f %f %d\n', N(i,1),N(i,2),k);
    end
end
fclose(fid);
toc;