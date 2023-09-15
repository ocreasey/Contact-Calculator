%Export Surface Objects Center of Homogeneous Mass

%Written by Olivia Creasey, PhD student
%This version written September 14-15, 2023
%Based off of previous "Surface Center Homogeneous Mass" XTensions 
%Previous versions written July 2023, edited September 2023

%
%<CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Export Surface Objects Center of Homogeneous Mass" icon="Matlab">
%          <Command>MatlabXT::XT_OAC_ContactCalculator_ExportCenterOfHomogenousMass(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces"> 
%          <Item name="Export Surface Objects Center of Homogeneous Mass" icon="Matlab">
%            <Command>MatlabXT::XT_OAC_ContactCalculator_ExportCenterOfHomogeneousMass(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>

%Description
%Iterates through all surface objects in the scene and exports XYZ
%center of homogeneous mass coordinates to a text file. This XTension only
%exports the center of homogenous mass coordinate for the first component
%within a given surface object - it will not export center of mass
%coordinates for additional components. If the surface object contains 
%disconnected components that have been unified, the XTension will export
%the aggregate center of mass coordinate for the entire surface object. 
%This is part of the suite of XTensions written for the Contact Calculator
%analysis but can be used outside of the Contact Calculator analysis.
%Written for use with Imaris 9.5.0
%
function XT_OAC_ContactCalculator_ExportCenterOfHomogeneousMass(aImarisApplicationID)
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

str=sprintf('Export center of homogeneous mass values for all surfaces?');
choice=questdlg(str, 'Center of homogeneous mass values export',...
    'Analyze','Cancel','Analyze');
%Handle Response
switch choice
    case 'Analyze'
    case 'Cancel'
        return;
end

% Exporting center of homogeneous mass statistics to a text file
vFilename = split(string(vImarisApplication.GetCurrentFileName),".");
filename = strcat(vFilename(1),'_','Center_Of_Homogeneous_Mass.txt');
fileID = fopen(filename,'w');

for d = 1:vNumberOfSurfaces
    vAllSurfacesStatistics = vSurfacesList{d}.GetStatistics;
    vAllSurfacesStatNames = cell(vAllSurfacesStatistics.mNames);
    vAllSurfacesStatValues = vAllSurfacesStatistics.mValues;
    vIndex=strmatch('Center of Homogeneous Mass',vAllSurfacesStatNames);
    
    vCenterHomogeneousMass = vAllSurfacesStatValues(vIndex,:);

    fprintf(fileID, '%s\r\n', vNamesList{d});
    fprintf(fileID, 'Center of Homogeneous Mass um  \t%f \t%f \t%f \r\n',vCenterHomogeneousMass);
    
    fprintf(fileID,'\r\n');

end
fclose(fileID);
end