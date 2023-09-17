%Surface Object Contact Locator

%Written by Olivia Creasey, PhD student
%This version written September 15-17, 2023
%Previous "Single Object Contact Calculator" versions written
%July-September 2023


%
%<CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Surface Object Contact Locator" icon="Matlab">
%          <Command>MatlabXT::XT_OAC_ContactCalculator_SurfaceObjectContactLocator(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces"> 
%          <Item name="Surface Object Contact Locator" icon="Matlab">
%            <Command>MatlabXT::XT_OAC_ContactCalculator_SurfaceObjectContactLocator(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%Description
%This is intended for use as part of the Contact Calculator analysis. It
%was specifically written to enable the user to find the surface areas and
%locations of cell-BM contacts in an image of a pancreatic islet, although
%it may have other uses in other applications of the Contact Calculator
%analysis. It also counts the discontinuous regions of contact that 
%particular cell has with a single source of BM (in this analysis, regions 
%of contact separated by a gap of more than ~1 um are considered to be 
%discontinous). 
%This XTension finds all of the surface areas of contact that a single
%surface object in the scene has with other tissue components, as 
%represented in the form of (cell) masks.  
%Start with a scene with surface objects and masks of all tissue components
%in a single channel as created using the Masks From Surfaces XTension. 
%For best performance, have only the surface object(s) that you want to
%analyze in this way, for speed/efficiency of analysis.
%A good way to do this is to create a new file in which you have deleted
%other surface objects OR create a new Imaris file with just the mask 
%channel and then use the surface object wizard to generate just the few 
%surface objects desired.
%For speed of calculation, have no other surface objects or channels in the
%file.
%Select as many surface objects in the scene as should be analyzed; the
%XTension will iterate through all of the selected surface objects and
%export individual text files for each surface object selected listing the
%identities of all of the contacts associated with that surface object, as
%well as the size, center of homogeneous mass (location) and number of
%discontinuous segments of each contact. 
%The user is also offered the option to generate unsmoothed or smoothed 
%surface objects that represent the contacts and that are retained in the
%Imaris file. 


function XT_OAC_ContactCalculator_SurfaceObjectContactLocator(aImarisApplicationID)
% connect to Imaris interface
if ~isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
    javaaddpath ImarisLib.jar
    vImarisLib = ImarisLib;
    if ischar(aImarisApplicationID)
        aImarisApplicationID = round(str2double(aImarisApplicationID));
    end
    vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
else
    vImarisApplication = aImarisApplicationID;
end

%We need the filename for the Imaris file so that we can write the
%computed results to a file in the same directory
vFilename = split(string(vImarisApplication.GetCurrentFileName),".");

% the user has to create a scene with some surfaces
vSurpassScene = vImarisApplication.GetSurpassScene;
if isequal(vSurpassScene, [])
    msgbox('Please create some Surfaces in the Surpass scene!');
    return;
end

vNumberOfChannels = vImarisApplication.GetDataSet.GetSizeC;
if vNumberOfChannels ~= 1
    msgbox('Please use a file that contains only one channel containing surface object masks');
    return;
end

vAnalysisModeList = {'Create unsmoothed contacts','Create smoothed contacts','Export data to text file'};
%Create Dialog box and allow user to choose the Primary Surfaces
[vPairMode, vOkMode] = listdlg('ListString',vAnalysisModeList,'SelectionMode','multiple',...
    'ListSize',[250 150],'Name','Surface Object Contact Locator','InitialValue',1, ...
    'PromptString',{'Please select all desired analysis outputs'});



if vOkMode == 0
    return
end
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
        vSurfacesList{vNumberOfSurfaces} = vImarisApplication.GetFactory.ToSurfaces(vDataItem);
        vNamesList{vNumberOfSurfaces} = char(vDataItem.GetName);
    end
end

if vNumberOfSurfaces<1
    msgbox('Please create at least 1 surface object!');
    return;
end


%%
%Choose the Primary Surfaces

%Create Dialog box and allow user to choose the Primary Surfaces
[vPair, vOk] = listdlg('ListString',vNamesList,'SelectionMode','multiple',...
    'ListSize',[250 150],'Name','Surface Object Contact Locator','InitialValue',1, ...
    'PromptString',{'Please select surfaces to analyze'});

if vOk == 0
    return
end


vPrimarySurface = cell(size(vPair));
vPrimarySurfaceList = cell(size(vPair));
for f = 1:length(vPair)
    vPrimarySurface{f} = vSurfacesList{vPair(f)};
    vPrimarySurfaceList{f} = char(vPrimarySurface{f}.GetName);
end
  

str=sprintf('Locate surface contacts for %s\n', vPrimarySurfaceList{:});
choice=questdlg(str, 'Locate surface contacts',...
    'Analyze','Cancel','Analyze');
%Handle Response
switch choice
    case 'Analyze'
    case 'Cancel'
        return;
end


%%
%Get Image Data parameters
vImarisDataSet0 = vImarisApplication.GetImage(0).Clone;
vImarisDataSet = vImarisApplication.GetImage(0).Clone;
vImarisDataSet1 = vImarisApplication.GetImage(0).Clone;


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

%Convert to 32bit
vFloatType = vImarisDataSet.GetType.eTypeFloat;
vImarisDataSet.SetType(vFloatType);
vImarisDataSet0.SetType(vFloatType);
vImarisDataSet1.SetType(vFloatType);

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

%%
vImarisApplication.SetVisible(0);
ip = vImarisApplication.GetImageProcessing;
Newsurfaces = vImarisApplication.GetFactory;
vRGBA=[255,255,255, 0];%for yellow
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]); % combine different components (four bytes) into one integer
%%
for PrimIter = 1:length(vPrimarySurface)
    
    vProgressDisplay = waitbar(0, 'Distance Transform: Preparation');

    % Create new channels for the analysis
    vImarisDataSet.SetSizeC(vNumberOfChannels + 2);
    
    

    %For a selected primary surface, mask the surface and then create a
    %distance transform
    %It seems that you have to mask first before you can do the distance
    %transform

    vMaskDataSetPrimary = vPrimarySurface{PrimIter}.GetMask( ...
            vDataMin(1), vDataMin(2), vDataMin(3), ...
            vDataMax(1), vDataMax(2), vDataMax(3), ...
            vDataSize(1), vDataSize(2), vDataSize(3), 0);
    
    %While looping in z, find out which z-slices actually contain the
    %primary surface, to save loops later while iterating through secondary
    %surfaces
    
    %Also determine which z and y slices actually contain the primary
    %surface, to reduce the calculations performed later
    
    
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
        
        if nnz(vSlice) == 0 && vPrimUpperZ == 0
            vPrimLowerZ = vIndexZ;           
        elseif nnz(vSlice) ~= 0
            vPrimUpperZ = vIndexZ;
        else
            continue
        end
        %This step puts the mask from Primary Surface into a channel, in
        %preparation for the distance transform

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
        
        if nnz(vSlice) == 0 && vPrimUpperY==0
            vPrimLowerY = vIndexY;
        elseif nnz(vSlice) ~= 0
            vPrimUpperY = vIndexY;
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
    
    vROI = [vStartX,vStartY,vStartZ,0,vEndX,vEndY,vEndZ,0];

    %create the distance transform for outside of the Primary Surface
    waitbar(0.5, vProgressDisplay, 'Distance Transform: Calculation');
    ip.DistanceTransformChannel(vImarisDataSet, vNumberOfChannels, 1, false); %puts it in channel 2
    waitbar(1, vProgressDisplay);
    close(vProgressDisplay);
    
    %Set the working dataset to vImarisDataSet only after adding the
    %distance transform channel to vImarisDataSet because (anecdotally) it 
    %is computationally faster to edit datasets that are not the working
    %dataset
    vImarisApplication.SetDataSet(vImarisDataSet);

    
    %Create a new folder object for new surfaces if requested by the user
    %Newsurfaces = vImarisApplication.GetFactory;
    if ismember(1,vPairMode)
        resultUnsmoothed = Newsurfaces.CreateDataContainer;
        resultUnsmoothed.SetName(sprintf('Unsmoothed Surfaces from %s',string(vPrimarySurfaceList{PrimIter})));
    end
    
    if ismember(2,vPairMode)
        resultSmoothed = Newsurfaces.CreateDataContainer;
        resultSmoothed.SetName(sprintf('Smoothed Surfaces from %s',string(vPrimarySurfaceList{PrimIter})));
    end    
    
    %Create a surface object containing the voxels at a distance of 1-1.5
    %voxel spacings from the outer surface of the Primary Surface
    vNewSurface = ip.DetectSurfacesWithUpperThreshold(vImarisDataSet,vROI,1,0,0,true,false,LowerThreshold,true,false,UpperThreshold,'');
    %This is a really useful method for identifying the voxels of interest
    %from a dataset stored as a channel in Imaris. If there are
    %discontinuous regions, they will be recognized as disconnected
    %components in the resulting surface object.
    
    vContactIntensityMatrix = zeros(vDataSize(1), vDataSize(2), vDataSize(3),'uint32');
    
    vProgressDisplay = waitbar(0, 'Generating Masks');
    
    %Create a matrix with zeros everywhere except in the voxels at a
    %distance of 1-1.5 voxel spacings from the outer surface of the primary
    %surface where the values are 1
    vMaskSurface = vNewSurface.GetMask( ...
            vDataMin(1), vDataMin(2), vDataMin(3), ...
            vDataMax(1), vDataMax(2), vDataMax(3), ...
            vDataSize(1), vDataSize(2), vDataSize(3), 0);
        
        
    for vIndexZ = vStartZ:vEndZ
        vSlice = vMaskSurface.GetDataSliceBytes(vIndexZ-1, 0, 0);
        if nnz(vSlice)==0
            continue %If the z-slice doesn't contain the Secondary Surface, just continue
        else
            %Create a matrix containing the intensity values from the 
            %surface objects mask channel
            vMaskIntensity = vImarisDataSet.GetDataSliceFloats(vIndexZ-1,0,0);
            vSlice = vSlice ==1;
            vSlice = uint32(vSlice).*uint32(vMaskIntensity);
            %Create a matrix with zeros everywhere except in the voxels at
            %a distance of 1-1.5 voxel spacings from the outer surface of
            %the primary surface were the values are the intensity values
            %from the surface objects mask channel
            vContactIntensityMatrix(:,:,vIndexZ) = vContactIntensityMatrix(:,:,vIndexZ)+uint32(vSlice);
            
        end
        waitbar((1+vIndexZ-vStartZ)/(vEndZ-vStartZ), vProgressDisplay);
    end
    close(vProgressDisplay);
 
    %Only some intensity values are used in the Masks channel, so we need to know what they are
    vContactIntensityVector = reshape(vContactIntensityMatrix,[],1);
    vChannelIntensityList = unique(vContactIntensityVector);
    
%%
    
    vProgressDisplay = waitbar(0, 'Producing Surface Objects');
    
    filenameprep = string(vPrimarySurfaceList{PrimIter});
    filename = strcat(vFilename(1),'_ContactLocator_',filenameprep,'.txt');
    fileID = fopen(filename,'w');
    fprintf(fileID, strcat(filenameprep, '\r\n'));
    
    for b = 2:length(vChannelIntensityList)
        %Start counting at 2 to ignore '0' voxels
        
        %For each contact, create a mask that will be used to generate
        %surfaces
        %First, a mask of '1's
        vSecondaryContactMatrix = vContactIntensityMatrix == vChannelIntensityList(b);
        %Next, a mask of '500's (this gives good results in our
        %application, but may need to be tailored for another user's
        %application)
        vSecondaryContactMatrixEstimated = vSecondaryContactMatrix*500;
  
        %This generates the ContactSurface which is the one we will use to
        %calcualte size of contact and location of contact
        %We already have vMaskSurface which is a DataSet, so we are going
        %to reuse it since we must turn vSecondaryContactMatrix into a
        %DataSet
        
        for vIndexZ = vStartZ:vEndZ
            %Put the secondary contact matrix into a DataSet
            vMaskSurface.SetDataSliceFloats(vSecondaryContactMatrix(:,:,vIndexZ),vIndexZ-1,0,0);   
            %Also put the 500 secondary contact matrix into a channel
            vImarisDataSet1.SetDataSliceFloats(vSecondaryContactMatrixEstimated(:,:,vIndexZ),vIndexZ-1,0,0);
        end

        
        %Create an empty surface, and then add a surface object using a
        %mask
        %This approach allows us to create a surface object in which
        %discontinuous surfaces are all "unified" into a single surface
        %object, avoiding any concerns about disconnected components
        vContactSurface = vImarisApplication.GetFactory.CreateSurfaces();
        vImarisApplication.GetFactory.ToSurfaces(vContactSurface);%To do this, mask values must be 1
        vContactSurface.AddSurface(vMaskSurface,0); 

        
        if ismember(1,vPairMode)
            vContactSurface.SetName(sprintf('%s Unsmoothed Surface %s',string(vPrimarySurfaceList{PrimIter}),string(vChannelIntensityList(b))));
            vContactSurface.SetColorRGBA(vRGBA);
            resultUnsmoothed.AddChild(vContactSurface, -1);
        end
        
        
        %This generates the SmoothedContactSurface which is the one we
        %will use to calculate the number of disconnected contacts (in our
        %test datasets with voxel size 0.1x0.1x0.1um, this set of 
        %parameters including the voxel intensity value of 500 detects 
        %contacts as disconnected if they are >~1um apart but not if they
        %are <1um apart
        vSmoothedContactSurface = ip.DetectSurfaces(vImarisDataSet1,vROI,0,2*MAXvoxelspacing,0,false,1,'');
        vContactSurfaceNumberOfContacts = vSmoothedContactSurface.GetNumberOfSurfaces();
        
        if ismember(2,vPairMode)
            vSmoothedContactSurface.SetName(sprintf('%s Smoothed Surface %s',string(vPrimarySurfaceList{PrimIter}),string(vChannelIntensityList(b))));
            vSmoothedContactSurface.SetColorRGBA(vRGBA);
            resultSmoothed.AddChild(vSmoothedContactSurface, -1);
        end  
        
        vStatistics = vContactSurface.GetStatistics;
        vStatisticsNames = cell(vStatistics.mNames);
        vStatisticsValues = vStatistics.mValues;
        vStatisticsIndex1=strmatch('Number of Voxels', vStatisticsNames);
        vNumberOfVoxels = vStatisticsValues(vStatisticsIndex1,:);
        vStatisticsIndex2=strmatch('Center of Homogeneous Mass X',vStatisticsNames);
        vCenterOfHomogeneousMassX = vStatisticsValues(vStatisticsIndex2(end),:);
        vStatisticsIndex3=strmatch('Center of Homogeneous Mass Y',vStatisticsNames);
        vCenterOfHomogeneousMassY = vStatisticsValues(vStatisticsIndex3(end),:);
        vStatisticsIndex4=strmatch('Center of Homogeneous Mass Z',vStatisticsNames);
        vCenterOfHomogeneousMassZ = vStatisticsValues(vStatisticsIndex4(end),:);

    
        fprintf(fileID,'\r\n');
        fprintf(fileID,'%s\r\n', string(vChannelIntensityList(b)));
        fprintf(fileID,'Number of Voxels  \t%f \r\n',vNumberOfVoxels);
        fprintf(fileID,'Center of Homogeneous Mass um  \t%f \t%f \t%f \r\n',vCenterOfHomogeneousMassX,vCenterOfHomogeneousMassY,vCenterOfHomogeneousMassZ);
        fprintf(fileID,'Number of Contacts \t%f \r\n',vContactSurfaceNumberOfContacts);
        
        waitbar(b/length(vChannelIntensityList), vProgressDisplay);
        
    end
    
    fclose(fileID);
    
    if ismember(1,vPairMode)
        vImarisApplication.GetSurpassScene.AddChild(resultUnsmoothed, -1);
    end
    
    if ismember(2,vPairMode)
        vImarisApplication.GetSurpassScene.AddChild(resultSmoothed, -1);
    end
    
    
    close(vProgressDisplay);
    vImarisApplication.SetDataSet(vImarisDataSet0);
    vImarisDataSet.Crop(0,vDataSize(1),0,vDataSize(2),0,vDataSize(3),0,vNumberOfChannels,0,aSizeT);
    vMaskSurface.Dispose();
    vMaskDataSetPrimary.Dispose();
    
end
    
vImarisApplication.SetDataSet(vImarisDataSet);    
vImarisApplication.SetVisible(1);    
vImarisDataSet.Dispose();
vImarisDataSet0.Dispose();
vImarisDataSet1.Dispose();

end