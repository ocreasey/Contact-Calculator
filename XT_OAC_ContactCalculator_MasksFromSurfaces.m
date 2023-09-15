%Export Masks from Surfaces

%Written by Olivia Creasey, PhD Student.
%This version written September 14, 2023
%Previous version written January 2023

%
%<CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Export Masks from Surfaces" icon="Matlab">
%          <Command>MatlabXT::XT_OAC_ContactCalculator_MasksFromSurfaces(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Export Masks from Surfaces" icon="Matlab">
%            <Command>MatlabXT::XT_OAC_ContactCalculator_MasksFromSurfaces(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%Description
%This XTension will create masks within a single channel for all 
%individual surface objects in the scene. The intensity value for each mask
%will correspond to the suface object's number identifier. The surface
%objects must be named something like "Surface 1" with a word and then a
%number separated by a space. The masks channel can be exported as a tiff
%from Imaris after the channel is generaged by this XTension.
%This XTension should only be used if the image has not been downsampled or
%upsampled since the surface objects were generated - the masks will be
%overlapping if more than one surface object occupies a given voxel.
%This is part of the suite of XTensions written for the Contact Calculator
%analysis but can be used outside of the Contact Calculator analysis.
%Written for use with Imaris 9.5.0

function XT_OAC_ContactCalculator_MasksFromSurfaces(aImarisApplicationID)
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

vImarisApplication.SetVisible(0);

%%
%Get Image Data parameters
vImarisDataSet = vImarisApplication.GetDataSet.Clone;

vNumberOfChannels = vImarisDataSet.GetSizeC;

vDataMin = [vImarisDataSet.GetExtendMinX, vImarisDataSet.GetExtendMinY, vImarisDataSet.GetExtendMinZ];
vDataMax = [vImarisDataSet.GetExtendMaxX, vImarisDataSet.GetExtendMaxY, vImarisDataSet.GetExtendMaxZ];
vDataSize = [vImarisDataSet.GetSizeX, vImarisDataSet.GetSizeY, vImarisDataSet.GetSizeZ];
aSizeT = vImarisApplication.GetDataSet.GetSizeT;

% Create a new matrix to hold the result 
vImage = zeros(vDataSize(1), vDataSize(2), vDataSize(3),'uint32');

vProgressDisplay = waitbar(0, 'Generating Masks');

for SurfaceIter = 1:vNumberOfSurfaces 
    
    vMaskSurface = vSurfacesList{SurfaceIter}.GetMask( ...
            vDataMin(1), vDataMin(2), vDataMin(3), ...
            vDataMax(1), vDataMax(2), vDataMax(3), ...
            vDataSize(1), vDataSize(2), vDataSize(3), 0);
    
    MaskIntensity = split(vNamesList{SurfaceIter},' ');
    vMaskIntensity = str2double(string(MaskIntensity(2)));
        
    for vIndexZ = 1:vDataSize(3)
        %Mask for Primary surface (inside=1)
        vSlice = vMaskSurface.GetDataSliceBytes(vIndexZ-1, 0, 0);
        vSlice = vSlice == 1;
        vSlice = vSlice*vMaskIntensity;
        vImage(:,:,vIndexZ) = vImage(:,:,vIndexZ)+uint32(vSlice);
    end
    waitbar(SurfaceIter/vNumberOfSurfaces, vProgressDisplay);
end

vImarisDataSet.SetSizeC(vNumberOfChannels+1);

for vIndexZ = 1:vDataSize(3)
    %Discrepancy in z index counting is because Imaris starts counting at 0
    %and Matlab starts counting at 1
    vImarisDataSet.SetDataSliceFloats(vImage(:,:,vIndexZ),vIndexZ-1,vNumberOfChannels,0);
end
vImarisApplication.SetDataSet(vImarisDataSet);
vImarisApplication.SetVisible(1);
close(vProgressDisplay);
msgbox('Masks for all surfaces generated!');

end
