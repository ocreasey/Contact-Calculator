%Surface Objects Distance from a Single Surface 

%Written by Olivia Creasey, PhD student
%This version written September 14, 2023
%Previous version written August 30, 2023

%
%<CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Export Surface Objects Distance from a Single Surface" icon="Matlab">
%          <Command>MatlabXT::XT_OAC_ContactCalculator_DistanceFromSurface(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Export Surface Objects Distance from a Single Surface" icon="Matlab">
%            <Command>MatlabXT::XT_OAC_ContactCalculator_DistanceFromSurface(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%Description
%This XTension allows the user to export to a text file the median, mean, 
%object center of mass, minimum, and maximum distances for all surface 
%objects in a scene from a single chosen surface object. These are all 
%statistics that are automatically calculated for a given surface object by
%Imaris once a distance transform has been generated from the single chosen
%surface object. 
%Surfaces at a minimum distance of less than sqrt(2)*1 voxel dimension from
%the chosen surface object can be considered to be touching the object
%The original intention within the context of the Contact Calculator 
%analysis of pancreatic islet structure is to allow the user to export data
%related to cell distance from peri-islet BM.
%This is part of the suite of XTensions written for the Contact Calculator
%analysis but can be used outside of the Contact Calculator analysis.
%Written for use with Imaris 9.5.0

function XT_OAC_ContactCalculator_DistanceFromSurface(aImarisApplicationID)
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

if vNumberOfSurfaces<2
    msgbox('Please create at least 2 surface objects!');
    return;
end

%Create Dialog box and allow user to choose the single Primary Surface
[vPair, vOk] = listdlg('ListString',vNamesList,'SelectionMode','single',...
    'ListSize',[250 150],'Name','Distance from a Single Surface','InitialValue',1, ...
    'PromptString',{'Please select one Single Surface'});

if vOk == 0
    return
end

str=sprintf('Export distance values from %s?',vNamesList{vPair});
choice=questdlg(str, 'Distance values export',...
    'Analyze','Cancel','Analyze');
%Handle Response
switch choice
    case 'Analyze'
    case 'Cancel'
        return;
end

vPrimarySurface = vSurfacesList{vPair};
vPrimarySurfaceName = vNamesList{vPair};

% Create new channels for the analysis
vImarisDataSet = vImarisApplication.GetImage(0).Clone;
vDataMin = [vImarisDataSet.GetExtendMinX, vImarisDataSet.GetExtendMinY, vImarisDataSet.GetExtendMinZ];
vDataMax = [vImarisDataSet.GetExtendMaxX, vImarisDataSet.GetExtendMaxY, vImarisDataSet.GetExtendMaxZ];
vDataSize = [vImarisDataSet.GetSizeX, vImarisDataSet.GetSizeY, vImarisDataSet.GetSizeZ];
aSizeT = vImarisApplication.GetDataSet.GetSizeT;

%Convert to 32bit
%This step is critical when calculating a distance transform - otherwise
%the distance transform calculated is only 8-bit
vFloatType = vImarisDataSet.GetType.eTypeFloat;
vImarisDataSet.SetType(vFloatType);


vImarisDataSet.SetSizeC(vNumberOfChannels + 1);
%Turn off the Imaris display to reduce time and CPU/RAM usage
vImarisApplication.SetVisible(0);
    
vProgressDisplay = waitbar(0, 'Distance Transform: Preparation');

%For a selected primary surface, mask the surface and then create a
    %distance transform

vMaskDataSetPrimary = vPrimarySurface.GetMask( ...
       vDataMin(1), vDataMin(2), vDataMin(3), ...
       vDataMax(1), vDataMax(2), vDataMax(3), ...
       vDataSize(1), vDataSize(2), vDataSize(3), 0);
   
   
 %You have to turn the mask (a DataSet) into an array so that you can set the mask as a channel  
 for vIndexZ = 1:vDataSize(3)
        %Mask for Primary surface (inside=1)
        vSlice = vMaskDataSetPrimary.GetDataSubVolumeAs1DArrayBytes(0,0,vIndexZ-1, 0, 0,vDataSize(1), vDataSize(2),1);
        
        %Set the intensity values for the new channel in Imaris to be the
        %surface object mask in vMaskDataSetPrimary - The seeming
        %discrepancy in the counting of the channels and in the iteration 
        %through slices in z is because Imaris starts counting at 0 and
        %Matlab starts counting at 1
        vImarisDataSet.SetDataSubVolumeAs1DArrayBytes(vSlice, ...
            0,0,vIndexZ-1,vNumberOfChannels,0,vDataSize(1),vDataSize(2),1);
        
        waitbar(vIndexZ/vDataSize(3), vProgressDisplay);
 end
 
%create the distance transform for outside of the Primary Surface
waitbar(0.5, vProgressDisplay, 'Distance Transform: Calculation');
%Replace the surface object mask with a distance transform
vImarisApplication.GetImageProcessing.DistanceTransformChannel( ...
    vImarisDataSet, vNumberOfChannels, 1, false); 
waitbar(1, vProgressDisplay);
close(vProgressDisplay);

vImarisApplication.SetDataSet(vImarisDataSet);

%Make a file for the Primary Surface and write information about the
%Primary Surface to the file
filenameprep = regexprep(vPrimarySurfaceName,' ','_');
filename = strcat(vFilename(1),'_Distance_from_',filenameprep,'.txt');
fileID = fopen(filename,'w');
fprintf(fileID, strcat(filenameprep, '\r\n'));

for i = 1:vNumberOfSurfaces
    %skip exporting anything for the Primary Surface
    if i == vPair
        continue
    end
    
    vStatistics = vSurfacesList{i}.GetStatistics;
    vStatisticsNames = cell(vStatistics.mNames);
    vStatisticsValues = vStatistics.mValues;
    vIndex=strmatch('Intensity Median',vStatisticsNames);
    vIndex1 = strmatch('Intensity Mean',vStatisticsNames);
    vIndex3 = strmatch('Intensity Center',vStatisticsNames);
    vIndex4 = strmatch('Intensity Min',vStatisticsNames);
    vIndex5 = strmatch('Intensity Max',vStatisticsNames);
    
    vChannelIntensityMedian = vStatisticsValues(vIndex,:);
    vChannelIntensityMean = vStatisticsValues(vIndex1,:);
    vChannelIntensityCenter = vStatisticsValues(vIndex3,:);
    vChannelIntensityMin = vStatisticsValues(vIndex4,:);
    vChannelIntensityMax = vStatisticsValues(vIndex5,:);
    
    fprintf(fileID,'\r\n');
    fprintf(fileID, '%s\r\n', vNamesList{i});
    fprintf(fileID,'Median Distance: \t%f\r\n', vChannelIntensityMedian(end));
    fprintf(fileID,'Mean Distance: \t%f\r\n', vChannelIntensityMean(end));
    fprintf(fileID,'Center Distance: \t%f\r\n', vChannelIntensityCenter(end));
    fprintf(fileID,'Min Distance: \t%f\r\n', vChannelIntensityMin(end));
    fprintf(fileID,'Max Distance: \t%f\r\n', vChannelIntensityMax(end));
    
end
fclose(fileID);
vImarisApplication.SetVisible(1);
 
end