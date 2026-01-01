
%Note: Run this code section by section to really understand the processing
%steps. There are several steps that are currently manually adjusted and
%require user iteration / understanding.

%This code is provided as is, as a part of a research project. 
%If you intend to use this code, review it thoroughly. Parts of the code
%are specific to the dataset provided in the sample files, hardcoded
%herein.

%Find the sample images/data in XXXXXXXXXXXXXXX.
%You will need the files:
%   BGsubtract.mat & Delta.mat (for SkipGradientIntegration=false)
%   Sinograms.mat (for SkipGradientIntegration=true)

%CC-BY Fernando Zigunov, 2026

%% 1. User Parameters
clear; clc; close all;

SkipGradientIntegration = false; %If true, skips Section (2) in this code (only for demo purposes)
 
%During a run, the background can shift as the components vibrate due to
%flow entrainment or thermal expansion of the light source backilluminating the BG.
%**Background subtraction requires no-refraction regions in the field of view in all images!**
%Select Background shift subtraction method:
BGsubtractionMethod = 'average'; %Simplest (use first). Just find average displacement by using all no-refraction vectors.
%BGsubtractionMethod = 'planefit'; %More complex (try only if 'average' fails; consider using if telecentric lens has low f-number; not used here). Fits a plane to the X and Y displacements. This will help correct for warping if the optics have barrel distortion.

%Experimental parameters
fps=20; %Camera FPS during acquisition
rpm=2; %Nozzle rotation RPM - used to determine pose angle for each pose

vectorSkip=12/4; %xcorr skips pixels which affects scaling
pxpermm = 9.9; %Scaling from calibration image
vecpermm=pxpermm/(vectorSkip); 

%Scaling factors 
DistanceSchlierenObjectToBG = 13.25*0.0254; %meters
GladstoneDale = 2.3e-4; %m^3/Kg
pxSize = 4.5e-6; %Camera pixel size in meters
Magnification = (891*pxSize/90e-3); %magnification factor (this example, a 90mm reference is 891px -> 9.9pxpermm)

ScalingFactor = pxSize/(DistanceSchlierenObjectToBG*Magnification*GladstoneDale); %This is the scaling constant from eq 6 in Zigunov 2025 Experiments in Fluids


%% 2. Integration of vector displacement fields into integrated gradients
% This step requires a gradient inversion solver (or a Poisson solver)
% Here we use our open-source OSMODI solver (https://github.com/3dfernando/pressure-osmosis)
% If you have trouble installing the CUDA version, consider using the CPU-only version (slower)

%To quick-start (skip this step), load the already-integrated fields in the file "Sinograms.mat"

if SkipGradientIntegration
    disp('Loading sample sinograms...')
    load('Sinograms.mat'); %File contains sinograms [1 variable 463x306x1200 in sample]
    disp('Complete!')
else
    disp('Loading sample vector displacements...')
    load('Delta.mat'); %File contains vector displacements [2 variables 534x334x1200 in sample data]

    disp('Processing vectors into sinograms...')
    skipImages = 1; %Change to do a quick processing skipping images to see the rotation more clearly
    for i=1:skipImages:size(DX,3)
        deltaX=DX(:,:,i); deltaY=DY(:,:,i);
        
        if i==1
            try
                load('BGsubtract.mat'); %Loads mask
            catch
                figure;
                imagesc(deltaY); title('Select Cropping Rectangle');
                R=getrect();
                R=round(R);
                deltaY_test=deltaY(R(2):(R(2)+R(4)),R(1):(R(1)+R(3)));
                
                imagesc(deltaY_test); 
                title('Select Area where Flow is Present Inside');
                daspect([1 1 1]);
                PBG=roipoly();
                close all;
    
                save('BGsubtract.mat','PBG','R');
            end
        end
        %Crops
        deltaX=deltaX(R(2):(R(2)+R(4)),R(1):(R(1)+R(3)));
        deltaY=deltaY(R(2):(R(2)+R(4)),R(1):(R(1)+R(3)));
    
        %removes 0 (zero = masked out)
        deltaX(deltaY==0 & ~PBG)=nan; deltaY(deltaY==0 & ~PBG)=nan;
        
        %Removes BG motion with plane fit on background
        %Select BGmethod (in section 1) to define behavior
    
        if strcmpi(BGsubtractionMethod,'planefit')
            ValidBG=~PBG & ~isnan(deltaX);
            [X,Y]=ndgrid(1:size(deltaX,1),1:size(deltaY,2));
            Matrix=[ones(sum(ValidBG(:)),1), X(ValidBG), Y(ValidBG)];
            FitX=Matrix\deltaX(ValidBG);
            FitY=Matrix\deltaY(ValidBG);
        
            deltaX_BG=FitX(1) + FitX(2) * X + FitX(3) * Y;
            deltaY_BG=FitY(1) + FitY(2) * X + FitY(3) * Y;
        
            deltaX=deltaX-deltaX_BG;
            deltaY=deltaY-deltaY_BG;
        elseif strcmpi(BGsubtractionMethod,'average')
            ValidBG=~PBG & ~isnan(deltaX);
            
            deltaX_BG=mean(deltaX(ValidBG),'omitnan');
            deltaY_BG=mean(deltaY(ValidBG),'omitnan');
              
            deltaX=deltaX-deltaX_BG;
            deltaY=deltaY-deltaY_BG;
        end
            
        %Fills missing vectors if vector field corrupted inside (from vector validation)
        deltaX=fillmissing(deltaX,'linear','EndValues','none');
        deltaY=fillmissing(deltaY,'linear','EndValues','none');
    
        % Solves for density
        [n,~]=OSMODI(deltaX,deltaY,zeros(size(deltaX))); %Solves with grid spacing = 1 (grid spacing not provided). Scaling will be applied later.
        n=n-mean(n(~PBG),'omitnan');
        %imagesc(n');daspect([1 1 1]); colorbar
    
        if i==1
            allSinograms=zeros(size(n,1),size(n,2),size(DX,3));
        end
    
        allSinograms(:,:,i)=n;
    
    
        imagesc(squeeze(allSinograms(:,:,i)')); colorbar; title(i); daspect([1 1 1])
        drawnow; pause(0.1)

        disp(['Image #' num2str(i) '/' num2str(size(DX,3))]);    
    end
    disp('Complete!')
end

%% 3. Find the center of rotation MANUALLY. 
%Run this loop manually and find the rows corresponding to the axis of
%rotation (yCenterStart, yCenterEnd). 
%Compensates for Rotation center tilted w.r.t. camera x-axis.

%===This needs to be improved!! ===
% In future: Measure axis rotation/wobble by taking a set of images with a
% physical, properly aligned reference shaft. 
yCenterStart=143; 
yCenterEnd=135;

widthStart=min(yCenterStart-1, size(allSinograms,2)-yCenterStart);
widthEnd=min(yCenterEnd-1, size(allSinograms,2)-yCenterEnd);

width = min(widthStart, widthEnd); %Easiest to do (inscribed cylinder); but if conical volume is of interest needs to be updated

skip=1; %Use to assess if sinogram is centered more quickly
for x=1:skip:size(allSinograms,1)
    yCenter = ((x-1)/(size(allSinograms,1)-1)) * (yCenterEnd - yCenterStart) + yCenterStart;
    newYvalues = (yCenter-width):(yCenter+width);

    if x==1
        allSinograms2=zeros(size(allSinograms,1),width*2+1,size(allSinograms,3));
    end        

    %Interpolates to get better centering and no jumps
    currentSinogram = squeeze(allSinograms(x,:,:));
    yVals = 1:size(allSinograms,2); fVals = 1:size(allSinograms,3);
    [yV, fV] = meshgrid(yVals,fVals);
    [yV2, fV2] = meshgrid(newYvalues,fVals);

    currentSinogram2 = interp2(yV,fV,currentSinogram',yV2,fV2);
    
    allSinograms2(x,:,:) = currentSinogram2';

%     %Shows the centerline to assess performance
%     imagesc(squeeze(allSinograms2(x,:,:))); title(x); hold on;
%     plot([0 size(allSinograms2,3)],[1 1]*width+0.5,'k--');        
%     drawnow; pause(0.01);

    disp(['Processing Sinogram ' num2str(x) '/' num2str(size(allSinograms,1))]);
end

%% 4. Performs filtered backprojection
totalFrames=size(allSinograms2,3);
nFullRotations=totalFrames/(fps*60/rpm);

theta=linspace(0,nFullRotations*360,totalFrames); %Builds theta array based on FPS and RPM

%Rearrange sinogram such that image is reconstructed at zero degrees
zRotationAngle = 5.1; %degrees
firstFrame = round((zRotationAngle / 360) * (totalFrames/nFullRotations));
newFrameOrder = [(firstFrame+1):totalFrames 1:firstFrame];

allSinograms2 = allSinograms2(:,:,newFrameOrder);

for x=1:skip:size(allSinograms,1)
    %Volume reconstruction slice-by-slice (parallel ray allows this)
    nSlice = -squeeze(allSinograms2(x,:,:)); %Negative sign required because vector displacement is inverted due to optics
    nSlice(isnan(nSlice))=0;

    outputSize = 2*floor(size(nSlice,1)/(2)); %Make the biggest volume allowed by iRadon  
    nVolumeSlice = iradon(nSlice,theta,outputSize);
    
    if x==1
        xVol = (1:size(nVolumeSlice,1))-(size(nVolumeSlice,1)/2);
        yVol = (1:size(nVolumeSlice,2))-(size(nVolumeSlice,2)/2);
        [xV, yV] = ndgrid(xVol,yVol);
        rVol = sqrt(xV.^2+yV.^2);
        volSliceMask = ones(size(nVolumeSlice));
        cutoffRadius = (size(nVolumeSlice,1)/2) * 0.95; %Removes a bit of the edge effects
        volSliceMask(rVol<cutoffRadius) = 0; %Mask to eliminate regions outside the projection domain
    end

    nVolumeSlice(volSliceMask==1) = nan; %Removes regions outside projection domain

    if x==1
        nVolume=zeros(size(nVolumeSlice,1),size(nVolumeSlice,2),size(allSinograms,1)); %Presizes
    end

    nVolume(:,:,x) = nVolumeSlice;

    if mod(x,10)==0
        imagesc(nVolumeSlice');daspect([1 1 1]); colorbar
        title(['Plane ' num2str(x) ]);
        drawnow; pause(0.05)
    end
end

%% 5. Rescales to physical units
xx=((1:size(nVolume,1))-round(size(nVolume,1)/2))/vecpermm;
yy=((1:size(nVolume,2))-round(size(nVolume,2)/2))/vecpermm;
zz=((1:size(nVolume,3))-round(size(nVolume,3)/2))/vecpermm;

rhoVolume=nVolume*ScalingFactor;
rhoInf=rhoVolume(35,129,450); %Constant shift at ambient (USER CHOSEN)
rhoVolume=rhoVolume-rhoInf + 1.2; %Adds density at ambient

lapRhoVolume = del2(rhoVolume,(xx(2)-xx(1))/1e3);
[dRho_dx,dRho_dy,dRho_dz] = gradient(rhoVolume,(xx(2)-xx(1))/1e3);
gradRhoMag = sqrt(dRho_dx.^2+dRho_dy.^2+dRho_dz.^2);

%Removes nans with zeros as the vtk doesn't like the nans
dRho_dx(isnan(dRho_dx))=0;
dRho_dy(isnan(dRho_dy))=0;
dRho_dz(isnan(dRho_dz))=0;
lapRhoVolume(isnan(lapRhoVolume))=0;
gradRhoMag(isnan(gradRhoMag))=0;
rhoVolume(isnan(rhoVolume))=0;

%Shows volume with Matlab's native volshow (just to see structure of flow)
volshow((rhoVolume-0.2).^8)

%% 6. Saves processed files and VTK for visualization
save('4-Nozzle Demo.mat','-v7.3');

%Exports VTK for Paraview visualization
vtkwrite('4-Nozzle Demo.vtk','RECTILINEAR_GRID',xx,yy,zz,'Scalars','rho[Kg_m3]',rhoVolume, ...
    'Scalars','LapRho[Kg_m5]',lapRhoVolume,'Scalars','|GradRho|[Kg_m4]',gradRhoMag);

