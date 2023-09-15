%Contact Calculator - Main Downsampled

%Written by Olivia Creasey, PhD Student.
%This version written September 14-15, 2023
%Previous versions written before November 2020, with additional edits 
%October 2021 and July 2023

%
%<CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Contact Calculator - Main (Downsampled)" icon="Matlab">
%          <Command>MatlabXT::XT_OAC_ContactCalculator_MainDownsampled(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Contact Calculator - Main (Downsampled)" icon="Matlab">
%            <Command>MatlabXT::XT_OAC_ContactCalculator_MainDownsampled(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%Description
%This XTension finds, measures, and exports measurements for surface areas 
%of contact between surface objects in an Imaris Scene. 
%This XTension was based on the algorithm in the Surface-Surface Contact
%Area XTension. The core principles of the algorithm remain, but the
%XTension no longer creates new surfaces corresponding to the calculated
%surface areas of contact between surface objects - it merely counts the
%voxels involved in those surface areas of contact and exports the data to
%text files. 
%This XTension expects and requires near-isotropic voxel dimensions.
%The updated algorithm automatically finds all of the surface objects 
%(Secondary Surfaces) that touch a user-specified surface object 
%(Primary Surface), determines whether or not the Primary Surface
%is wholly embedded in the image, counts the voxels that are involved in
%the surface area of contact between the Primary Surface and each Secondary
%Surface (accounting for voxels counted multiple times due to
%down-sampling/digital calculations), and exports data about the morphology
%of the Primary Surface and the sizes of all of the contacts with Secondary
%Surfaces to a set of individual text files, one per Primary Surface. 
%This version of Contact Calculator is suitable for use if the image 
%has been downsampled after the generation of the surface objects. In
%general, we do not recommend changing the image sampling after generation
%of the surface objects. This version will also work for images that have
%not been downsampled, but it will be less efficient in its calculations.
%The XTension automatically iterates through the full list of
%user-specified Primary Surfaces. 
%The XTension has been written to perform calculations as quickly and
%efficiently as we could, requiring a local minimum of computational power 
%and RAM.
%Written for use with Imaris 9.5.0

function XT_OAC_ContactCalculator_MainDownsampled(aImarisApplicationID)
% Set up connection	between	Imaris and MATLAB
if	isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
    vImarisApplication = aImarisApplicationID;
else
    % connect to Imaris interface
    javaaddpath ImarisLib.jar
    vImarisLib = ImarisLib;
    if	ischar(aImarisApplicationID)
        aImarisApplicationID = round(str2double(aImarisApplicationID));
    end
    vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
end

% the user has to create a scene with some surfaces
vSurpassScene = vImarisApplication.GetSurpassScene;
if isequal(vSurpassScene, [])
    msgbox('Please create some Surfaces in the Surpass scene!');
    return;
end

%We need the filename for the Imaris file so that we can write the
%computed results to a file in the same directory
vFilename = split(string(vImarisApplication.GetCurrentFileName),".");
%%
% get all Surpass surfaces names
vSurfaces = vImarisApplication.GetFactory.ToSurfaces(vImarisApplication.GetSurpassSelection);
vSurfacesSelected = vImarisApplication.GetFactory.IsSurfaces(vSurfaces);

if vSurfacesSelected
    vScene = vSurfaces.GetParent;
else
    vScene = vImarisApplication.GetSurpassScene;
end
vNumberOfSurfaces = 0;
for vChildIndex = 1:vScene.GetNumberOfChildren
    vDataItem = vScene.GetChild(vChildIndex - 1);
    if vImarisApplication.GetFactory.IsSurfaces(vDataItem)
        vNumberOfSurfaces = vNumberOfSurfaces+1;
        %You have to cast these surfaces as surfaces even though they were 
        %just identified that were just identified as surfaces. I don't
        %know why, but you do
        vSurfacesList{vNumberOfSurfaces} = vImarisApplication.GetFactory.ToSurfaces(vDataItem);
        vNamesList{vNumberOfSurfaces} = char(vDataItem.GetName);
    end
end

if vNumberOfSurfaces<2
    msgbox('Please create at least 2 surfaces objects!');
    return;
end


%%
%Choose the Primary Surfaces

%Create Dialog box and allow user to choose the Primary Surfaces
[vPair, vOk] = listdlg('ListString',vNamesList,'SelectionMode','multiple',...
    'ListSize',[250 150],'Name','Surface-Surface Contact Area','InitialValue',1, ...
    'PromptString',{'Please select Primary Surfaces'});

if vOk == 0
    return
end


vPrimarySurface = cell(size(vPair));
vPrimarySurfaceList = cell(size(vPair));
for f = 1:length(vPair)
    vPrimarySurface{f} = vSurfacesList{vPair(f)};
    vPrimarySurfaceList{f} = char(vPrimarySurface{f}.GetName);
end
  

str=sprintf('Calculate surface contacts for %s\n', vPrimarySurfaceList{:});
choice=questdlg(str, 'Surface-Surface contact area analysis',...
    'Analyze','Cancel','Analyze');
%Handle Response
switch choice
    case 'Analyze'
    case 'Cancel'
        return;
end

%%
%Get Image Data parameters
vImarisDataSet0 = vImarisApplication.GetDataSet.Clone;
vImarisDataSet = vImarisApplication.GetDataSet.Clone;
vNumberOfChannels = vImarisDataSet.GetSizeC;

vDataMin = [vImarisDataSet.GetExtendMinX, vImarisDataSet.GetExtendMinY, vImarisDataSet.GetExtendMinZ];
vDataMax = [vImarisDataSet.GetExtendMaxX, vImarisDataSet.GetExtendMaxY, vImarisDataSet.GetExtendMaxZ];
vDataSize = [vImarisDataSet.GetSizeX, vImarisDataSet.GetSizeY, vImarisDataSet.GetSizeZ];
aSizeT = vImarisApplication.GetDataSet.GetSizeT;
Xvoxelspacing = (vDataMax(1)-vDataMin(1))/vDataSize(1);
Yvoxelspacing = (vDataMax(2)-vDataMin(2))/vDataSize(2);
Zvoxelspacing = (vDataMax(3)-vDataMin(3))/vDataSize(3);
%Account for the inherent anisotropy in voxel spacing in xy vs z
MAXvoxelspacing = max(Xvoxelspacing,Zvoxelspacing);
MINvoxelspacing = min(Xvoxelspacing,Zvoxelspacing);

%Set the thresholds to describe distances from the Primary Surface that
%correspond to one voxel
LowerThreshold=0.5*MINvoxelspacing;
UpperThreshold=1.5*MAXvoxelspacing;

ZLimit=((3*Xvoxelspacing)+100*Xvoxelspacing)/100;%test percent Xvoxelsize

%Convert to 32bit - necessary for proper distance transform calculation
vFloatType = vImarisDataSet.GetType.eTypeFloat;
vImarisDataSet.SetType(vFloatType);

%Convert data to isotropic voxels
if Zvoxelspacing>ZLimit
    endZadj=round(double(Zvoxelspacing/Xvoxelspacing*vDataSize(3)));
    str=sprintf('Warning this volume is does NOT have isotropic voxels.\n\n To properly run this application go to:\n Edit--""Resample3D"" and adjust Z-size to %d steps',endZadj);
    choice2=questdlg(str,'WARNING','Cancel','Cancel');
    %Handle Response
    switch choice2
        case 'Cancel'
            return;
    end
end


for PrimIter = 1:length(vPrimarySurface)

    %Use a distance transform to find the secondary surfaces that are 1
    %voxel distance from the Primary Surface
    vProgressDisplay = waitbar(0, 'Distance Transform: Preparation');

    % Create a new channel where the result will be sent
    vImarisDataSet.SetSizeC(vNumberOfChannels + 1);
    vImarisApplication.SetVisible(0);

    %For a selected primary surface, mask the surface and then create a
    %distance transform

    vMaskDataSetPrimary = vPrimarySurface{PrimIter}.GetMask( ...
            vDataMin(1), vDataMin(2), vDataMin(3), ...
            vDataMax(1), vDataMax(2), vDataMax(3), ...
            vDataSize(1), vDataSize(2), vDataSize(3), 0);
    
    %While looping in z, find out which z-slices actually contain the
    %primary surface, to save loops later while iterating through secondary
    %surfaces
    
    %Also determine which z and y slices actually contain the primary
    %surface, to reduce the calculations performed later
    
    %Also find if the primary surface is touching the edge of the image 
    vDistImageEdge = 1;
    vPrimLowerZ = 0;
    vPrimUpperZ = 0;
    vPrimLowerX = 0;
    vPrimUpperX = 0;
    vPrimLowerY = 0;
    vPrimUpperY = 0;

    for vIndexZ = 1:vDataSize(3)
        %Mask for Primary surface (inside=1)
        vSlice = vMaskDataSetPrimary.GetDataSubVolumeAs1DArrayBytes(0,0,vIndexZ-1, 0, 0,vDataSize(1), vDataSize(2),1);
        vSlice = vSlice == 1;
        if vIndexZ==1 || vIndexZ==vDataSize(3)
            if nnz(vSlice) ~= 0 %see if primary surface is touching the edge of the image in z
                vDistImageEdge = 0; 
            end
        end
        
        if nnz(vSlice) == 0 && vPrimUpperZ == 0
            vPrimLowerZ = vIndexZ;           
        elseif nnz(vSlice) ~= 0
            vPrimUpperZ = vIndexZ;
        else
            continue
        end

    end
    waitbar(1/6, vProgressDisplay);
    %Define the starting and ending z-slices for the loop through the
    %Secondary Surfaces, with a 3-voxel buffer on each side (or running
    %into the edges of the data set)
    vStartZ = max(1,vPrimLowerZ-3);
    vEndZ = min(vDataSize(3), vPrimUpperZ+3);
    vSizeZ = vEndZ - vStartZ;
    
    %Find edges in X
    for vIndexX = 1:vDataSize(1)
        vSlice = vMaskDataSetPrimary.GetDataSubVolumeAs1DArrayBytes(vIndexX-1, 0, vStartZ, 0, 0, 1, vDataSize(2), vSizeZ);
        vSlice = vSlice == 1;
        if vIndexX==1 || vIndexX==vDataSize(1)
            if nnz(vSlice) ~=0   %see if primary surface is touching the edge of the image in X
                vDistImageEdge=0;
            end
        end
        
        if nnz(vSlice) == 0 && vPrimUpperX==0
            vPrimLowerX = vIndexX;
        elseif nnz(vSlice) ~= 0
            vPrimUpperX = vIndexX;
        else
            continue
        end
    end
    waitbar(2/6, vProgressDisplay);
    
    vStartX = max(1,vPrimLowerX-3);
    vEndX = min(vDataSize(1), vPrimUpperX+3);
    vSizeX = vEndX - vStartX;
    
    
    %Find edges in Y
    for vIndexY = 1:vDataSize(2)
        vSlice = vMaskDataSetPrimary.GetDataSubVolumeAs1DArrayBytes(vStartX, vIndexY-1, vStartZ, 0, 0, vSizeX, 1, vSizeZ);
        vSlice = vSlice == 1;
        if vIndexY==1 || vIndexY==vDataSize(2)
            if nnz(vSlice) ~=0   %see if primary surface is touching the edge of the image in Y
                vDistImageEdge=0;
            end
        end
        
        if nnz(vSlice) == 0 && vPrimUpperY==0
            vPrimLowerY = vIndexY;
        elseif nnz(vSlice) ~= 0
            vPrimUpperY = vIndexY;
            %Set the intensity values in the channel only once the xyz ROI
            %has been identified, to increase efficiency
            vImarisDataSet.SetDataSubVolumeAs1DArrayBytes(vSlice, ...
            vStartX,vIndexY-1,vStartZ,vNumberOfChannels,0,vSizeX,1,vSizeZ);
        else
            continue
        end
        waitbar((2+(vIndexY+1)/vDataSize(2))/6, vProgressDisplay);
    end
    
    vStartY = max(1,vPrimLowerY-3);
    vEndY = min(vDataSize(2), vPrimUpperY+3);
    vSizeY = vEndY - vStartY;

    waitbar(0.5, vProgressDisplay, 'Distance Transform: Calculation');
    vImarisApplication.GetImageProcessing.DistanceTransformChannel( ...
        vImarisDataSet, vNumberOfChannels, 1, false);
    waitbar(1, vProgressDisplay);
    close(vProgressDisplay);

    vImarisApplication.SetDataSet(vImarisDataSet);

    vProgressDisplay = waitbar(0, 'Selecting Secondary Surfaces');

    %Identify the Secondary Surfaces physically in contact with the Primary
    %Surface
    vSecondarySurface={};
    for h=1:vNumberOfSurfaces
        if vPrimarySurface{PrimIter} == vSurfacesList{h}
            continue
        else
            vStatistics = vSurfacesList{h}.GetStatistics;
            vStatisticsNames = cell(vStatistics.mNames);
            vStatisticsValues = vStatistics.mValues;
            vStatisticsIndex1=strmatch('Intensity Min', vStatisticsNames);
            vMinDistanceValue = vStatisticsValues(vStatisticsIndex1(end));

            if vMinDistanceValue <= MAXvoxelspacing
                vSecondarySurface{end+1} = vSurfacesList{h};
            else
                continue
            end
            
        end
        waitbar(h/vNumberOfSurfaces, vProgressDisplay);
    end
    
    g=numel(vSecondarySurface);
    close(vProgressDisplay);

    %Get stats from the original primary surface
    vAllPrimaryStatistics = vPrimarySurface{PrimIter}.GetStatistics;
    vPrimaryStatNames = cell(vAllPrimaryStatistics.mNames);
    vPrimaryStatValues = vAllPrimaryStatistics.mValues;
    vIndex1=strmatch('Area', vPrimaryStatNames);
    vIndex2=strmatch('Volume', vPrimaryStatNames);
    vIndex3=strmatch('Sphericity', vPrimaryStatNames);
    vIndex4=strmatch('BoundingBoxOO', vPrimaryStatNames);
    vIndex5=strmatch('Ellipticity (prolate)',vPrimaryStatNames);
    vIndex6=strmatch('Ellipticity (oblate)',vPrimaryStatNames);
    vIndex7=strmatch('Center of Homogeneous Mass',vPrimaryStatNames);
    vPrimaryArea = vPrimaryStatValues(vIndex1,:);
    vPrimaryVolume = vPrimaryStatValues(vIndex2,:);
    vPrimarySphericity = vPrimaryStatValues(vIndex3,:);
    vPrimaryBoundingBox = vPrimaryStatValues(vIndex4,:);
    vPrimaryEllipticityP = vPrimaryStatValues(vIndex5,:);
    vPrimaryEllipticityO = vPrimaryStatValues(vIndex6,:);
    vCenterHomogeneousMass = vPrimaryStatValues(vIndex7,:);
    
    %That's all of the primary surface information sorted. Now for the
    %secondary surfaces
    

    %This is where we keep track of which voxels have been counted as part
    %of a contact, so that later we can identify overlaps
    %This is probably not necessary if there has been no
    %downsampling/upsampling 
    vXYSize = (vSizeX)*(vSizeY); 
    vROI_Cell = cell(vXYSize, vSizeZ+1);
    
    vResultValues=zeros(g); %This is where we save the voxel counts for each contact
    vProgressDisplay = waitbar(0, 'Surface Contact Area Calculation');
    %g is the total number of secondary surfaces
    for q=1:g
        vSecondaryMask = vSecondarySurface{q}.GetMask(vDataMin(1), vDataMin(2), vDataMin(3), vDataMax(1), vDataMax(2), vDataMax(3), vDataSize(1), vDataSize(2), vDataSize(3), 0);       
        
        for vIndexZ = vStartZ:vEndZ
            %Mask for secondary Surface (inside=1)
            vSlice2=vSecondaryMask.GetDataSubVolumeAs1DArrayFloats(vStartX,vStartY,vIndexZ-1,0,0,vSizeX,vSizeY,1);
            if nnz(vSlice2)==0
                continue %If the z-slice doesn't contain the Secondary Surface, just continue
            else
                vROI_index = vIndexZ-vStartZ+1;
                ch1=double(vImarisDataSet.GetDataSubVolumeAs1DArrayFloats(vStartX, vStartY, vIndexZ-1, vNumberOfChannels, 0, vSizeX, vSizeY, 1));
                %ch1 is from the distance transform channel for outside the
                %primary surface
                vSlice2 = vSlice2 == 1;% for inside secondary surface only
                ch2 = double(vSlice2);
                vData=(ch1.*ch2);
                for v = 1:vXYSize
                    if (LowerThreshold < vData(v)) && (vData(v) < UpperThreshold)%If a specific voxel is inside the Secondary Surface and directly adjacent to the Primary Surface
                        vROI_Cell{v,vROI_index} = [vROI_Cell{v,vROI_index},q];%keep track of contact voxel position, to adjust for overcounting
                        vResultValues(q)=vResultValues(q)+1; %Count the number of voxels for the contact
                    end
                end
            end
            
        end
        waitbar(q/g, vProgressDisplay);
        vSecondaryMask.Dispose();

    end
    
    %Find shared voxels and prepare to split evenly between the contacts that share
    %each voxel - will be subtracting vShared from vResultValues
    vShared = zeros(g);
    for s = 1:vSizeZ+1
        for t = 1:vXYSize
            if length(vROI_Cell{t,s})>1
                share = 1 - 1/length(vROI_Cell{t,s});
                for u = 1:length(vROI_Cell{t,s})
                    vShared(vROI_Cell{t,s}(u)) = vShared(vROI_Cell{t,s}(u)) + share;
                end
            end
        end
    end
    vResultValues = vResultValues - vShared;

         
    
    %Compute total calculated surface area
    vCalculatedSurfaceVoxels = sum(vResultValues,'all');
    
    %Make a file for the Primary Surface and write information about the
    %Primary Surface to the file
    filenameprep = regexprep(vPrimarySurfaceList{PrimIter},' ','_');
    filename = strcat(vFilename(1),'_ContactCalculator_',filenameprep,'.txt');
    fileID = fopen(filename,'w');
    fprintf(fileID, strcat(filenameprep, '\r\n'));
    fprintf(fileID, 'Surface Area um^2 \t%f\r\n',vPrimaryArea);
    fprintf(fileID, 'Volume um^3  \t%f \r\n',vPrimaryVolume);
    fprintf(fileID, 'Center of Homogeneous Mass um  \t%f \t%f \t%f \r\n',vCenterHomogeneousMass);
    fprintf(fileID, 'Sphericity   \t%f \r\n',vPrimarySphericity);
    fprintf(fileID, 'Bounding Box um  \t%f \t%f \t%f \r\n',vPrimaryBoundingBox);
    fprintf(fileID, 'Ellipticity (prolate) \t%f \r\n',vPrimaryEllipticityP);
    fprintf(fileID, 'Ellipticity (oblate) \t%f \r\n',vPrimaryEllipticityO);
    fprintf(fileID, 'Calculated Surface Voxels \t%f \r\n',vCalculatedSurfaceVoxels);
    fprintf(fileID, 'Distance to Border (Calculated): \t%f \r\n',vDistImageEdge);
    fprintf(fileID, 'Calculated Contact Area Voxels \tVoxel Size um^2 \t%f \r\n',(nthroot((Xvoxelspacing*Yvoxelspacing*Zvoxelspacing),3))^2);

    %Write the values for calculated contact area to the file
    for r = 1:g
        fprintf(fileID,'\r\n');
        fprintf(fileID, strcat(regexprep(string(char(vSecondarySurface{r}.GetName)),' ','_'),'\r\n'));
        fprintf(fileID,'%f\r\n',vResultValues(r));
    end          
    
    fclose(fileID);
    close(vProgressDisplay);
    vImarisApplication.SetDataSet(vImarisDataSet0);
    vImarisDataSet.Crop(0,vDataSize(1),0,vDataSize(2),0,vDataSize(3),0,vNumberOfChannels,0,aSizeT);
    vMaskDataSetPrimary.Dispose();
    clear vResultValues vSurfaceContact vAllSurfaceArea vSecondarySurface vSecondaryMask vMaskDataSetPrimary vSlice vSlice2 vData ch1 ch2 vROI_Cell vShared
    clear vPrimaryArea vPrimaryVolume vPrimarySphericity vPrimaryBoundingBox vPrimaryEllipticityP vPrimaryEllipticityO vCalculatedSurfaceVoxels vDistImageEdge
    clear vAllPrimaryStatistics vPrimaryStatNames vPrimaryStatValues vIndex1 vIndex2 vIndex3 vIndex4 vIndex5 vIndex6 vStatistics vStatisticsValues vStatisticsNames vStatisticsIndex1
    clear filenameprep filename fileID r s t u q h g vProgressDisplay vIndexX vIndexY vIndexZ vMinDistanceValue vXYSize
    clear vPrimLowerZ vPrimUpperZ vStartZ vEndZ vSizeZ vPrimLowerX vPrimUpperX vStartX vEndX vSizeX vPrimLowerY vPrimUpperY vStartY vEndY vSizeY vROI_Index share
end


vImarisApplication.SetVisible(1);
vImarisDataSet.Dispose();
vImarisDataSet0.Dispose();
end

