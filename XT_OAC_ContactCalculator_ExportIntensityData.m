%Export Surface Objects Intensity Data

%Written by Olivia Creasey, PhD student
%This version written September 14, 2023
%Based off of previous "Surface Intensity Values" XTensions 
%Previous versions written January 2020, edited March & November 2020


%
%<CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Export Surface Objects Intensity Data" icon="Matlab">
%          <Command>MatlabXT::XT_OAC_ContactCalculator_ExportIntensityData(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Export Surface Objects Intensity Data" icon="Matlab">
%            <Command>MatlabXT::XT_OAC_ContactCalculator_ExportIntensityData(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%

%Description
%Iterates through all surface objects in the scene and exports channel
%intensity data for all channels to a text file. The intensity data 
%exported are the data automatically generated for each surface object - 
%absolute maximum and minimum intensity values for each channel at any 
%voxel within the surface object, median and mean intensity values for each
%channel from among all of the voxels within the surface object, standard
%deviation of intensity values for each channel among all of the voxels
%within the surface object, and the sum of all intensity values for all
%voxels inside the surface object. A user can easily modify the XTension to
%omit export of certain statistics or to add the export of other
%statistics. 
%This is part of the suite of XTensions written for the Contact Calculator
%analysis but can be used outside of the Contact Calculator analysis.
%Written for use with Imaris 9.5.0

function XT_OAC_ContactCalculator_ExportIntensityData(aImarisApplicationID)

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

% the user has to create a scene with some surfaces
vSurpassScene = vImarisApplication.GetSurpassScene;
if isequal(vSurpassScene, [])
    msgbox('Please create some Surfaces in the Surpass scene!');
    return;
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
SurfacesList{vScene.GetNumberOfChildren} = [];
vNamesList{vScene.GetNumberOfChildren} = [];
for vChildIndex = 1:vScene.GetNumberOfChildren
    vDataItem = vScene.GetChild(vChildIndex - 1);
    if vImarisApplication.GetFactory.IsSurfaces(vDataItem)
        vNumberOfSurfaces = vNumberOfSurfaces+1;
        vSurfacesList{vNumberOfSurfaces} = vImarisApplication.GetFactory.ToSurfaces(vDataItem);
        vNamesList{vNumberOfSurfaces} = char(vDataItem.GetName);
    end
end

if vNumberOfSurfaces<1
    msgbox('Please create at least 1 surfaces object!');
    return;
end

vNamesList = vNamesList(1:vNumberOfSurfaces);
vSurfacesList = vSurfacesList(1:vNumberOfSurfaces);

str=sprintf('Export channel intensity values for all surfaces?');
choice=questdlg(str, 'Surface channel intensity values export',...
    'Analyze','Cancel','Analyze');
%Handle Response
switch choice
    case 'Analyze'
    case 'Cancel'
        return;
end

%Get Image Data parameters
vImarisDataSet = vImarisApplication.GetDataSet.Clone;
vNumberOfChannels = vImarisDataSet.GetSizeC;

% Exporting a bunch of information about the surfaces in the scene - the
% intensities of fluorescence in the various channels

vFilename = split(string(vImarisApplication.GetCurrentFileName),".");
filename = strcat(vFilename(1),'_','Channel_Intensity_Data.txt');
fileID = fopen(filename,'w');

for d = 1:vNumberOfSurfaces
    vAllSurfacesStatistics = vSurfacesList{d}.GetStatistics;
    vAllSurfacesStatNames = cell(vAllSurfacesStatistics.mNames);
    vAllSurfacesStatValues = vAllSurfacesStatistics.mValues;
    vIndex=strmatch('Intensity Median',vAllSurfacesStatNames);
    vIndex1 = strmatch('Intensity Mean',vAllSurfacesStatNames);
    vIndex3 = strmatch('Intensity Sum',vAllSurfacesStatNames);
    vIndex4 = strmatch('Intensity Min',vAllSurfacesStatNames);
    vIndex5 = strmatch('Intensity Max',vAllSurfacesStatNames);
    vIndex6 = strmatch('Intensity StdDev',vAllSurfacesStatNames);
    
    vChannelIntensityMedian = vAllSurfacesStatValues(vIndex,:);
    vChannelIntensityMean = vAllSurfacesStatValues(vIndex1,:);
    vChannelIntensitySum = vAllSurfacesStatValues(vIndex3,:);
    vChannelIntensityMin = vAllSurfacesStatValues(vIndex4,:);
    vChannelIntensityMax = vAllSurfacesStatValues(vIndex5,:);
    vChannelIntensityStd = vAllSurfacesStatValues(vIndex6,:);

    fprintf(fileID, '%s\r\n\r\n', vNamesList{d});
    for ch = 1:(vNumberOfChannels)
        fprintf(fileID,'Median Intensity Channel %d: \t%f\r\n', round(ch),vChannelIntensityMedian(ch));
    end
    for ch = 1:(vNumberOfChannels)
        fprintf(fileID,'Mean Intensity Channel %d: \t%f\r\n', round(ch),vChannelIntensityMean(ch));
    end
    for ch = 1:(vNumberOfChannels)
        fprintf(fileID,'Intensity Sum Channel %d: \t%f\r\n', round(ch),vChannelIntensitySum(ch));
    end
    for ch = 1:(vNumberOfChannels)
        fprintf(fileID,'Min Intensity Channel %d: \t%f\r\n', round(ch),vChannelIntensityMin(ch));
    end
    for ch = 1:(vNumberOfChannels)
        fprintf(fileID,'Max Intensity Channel %d: \t%f\r\n', round(ch),vChannelIntensityMax(ch));
    end
    for ch = 1:(vNumberOfChannels)
        fprintf(fileID,'Intensity StdDev Channel %d: \t%f\r\n', round(ch),vChannelIntensityStd(ch));
    end
    fprintf(fileID,'\r\n');

end
fclose(fileID);

end