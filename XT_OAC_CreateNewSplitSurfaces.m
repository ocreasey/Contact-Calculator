%Create New Split Surfaces

%Written by Olivia Creasey, PhD Student
%This version finalized September 14, 2023
%Original version written September 2019.  
%Originally written as a tool for the Contact Calculator analysis, but
%ultimately not used as part of the analysis pipeline
%
%<CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Create New Split Surfaces" icon="Matlab">
%          <Command>MatlabXT::XT_OAC_CreateNewSplitSurfaces(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Create New Split Surfaces" icon="Matlab">
%            <Command>MatlabXT::XT_OAC_CreateNewSplitSurfaces(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%Description
%This XTension is intended to help split specified connected surface objects after 
%they are generated through the Create Surfaces from Multicut XTension.
%This XTension will be particularly useful for splitting surfaces that are
%accidentally merged by the Multicut workflow and that are connected by
%only a few voxels. This XTension requires some setup, however. 
%The user must first create a mask of the surface that
%needs to be split and then use the Surfaces wizard to create surfaces off
%of that mask. The user should select no smoothing and use the Background
%Subtraction thresholding. Then enable Split Touching Objects with the
%appropriate Seed Points Diameter. At the end, there should be two new
%surfaces that are separated. In general, the volumes will be a little
%increased over the volume of the original surface. These new surfaces
%should then be individually duplicated to new, separate surface objects.
%Then, this XTension can be used to create new surfaces that are now split
%and approximately the same volume/shape as the original surface object. 
%This XTension was written for the Contact Calculator analysis but
%ultimately is not an integral or recommended component of the analysis 
%and therefore was not developed or used after 2019.
%Nevertheless, it may be useful in some way to others.
%Written for Imaris 9.5.0



function XT_OAC_CreateNewSplitSurfaces(aImarisApplicationID)

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
vSurfacesList{vScene.GetNumberOfChildren} = [];
vNamesList{vScene.GetNumberOfChildren} = [];
for vChildIndex = 1:vScene.GetNumberOfChildren
    vDataItem = vScene.GetChild(vChildIndex - 1);
    if vImarisApplication.GetFactory.IsSurfaces(vDataItem)
        vNumberOfSurfaces = vNumberOfSurfaces+1;
        vSurfacesList{vNumberOfSurfaces} = vImarisApplication.GetFactory.ToSurfaces(vDataItem);
        vNamesList{vNumberOfSurfaces} = char(vDataItem.GetName);
    end
end

if vNumberOfSurfaces<2
    msgbox('Please create at least 2 surfaces objects!');
    return;
end

vNamesList = vNamesList(1:vNumberOfSurfaces);
%%
%Choose the surfaces
%Choose how many surfaces to colocalize
vPrimarySurface=1;
vSecondarySurface=[];

%Create Dialog box and allow user to choose the Reference Position
vPair1 = [];
[vPair1, vOk] = listdlg('ListString',vNamesList,'SelectionMode','multiple',...
    'ListSize',[250 150],'Name','Create New Split Surface','InitialValue',[1,1], ...
    'PromptString',{'Please select Original Surface'});

if vOk == 0
    return
end

vPrimarySurface = vSurfacesList{vPair1(1)};
%Create Dialog box and allow user to choose the Reference Position
vPair2 = [];
[vPair2, vOk] = listdlg('ListString',vNamesList,'SelectionMode','multiple',...
    'ListSize',[250 150],'Name','Create New Split Surface','InitialValue',[1,1], ...
    'PromptString',{'Please select Wizard-Created Split Surfaces'});

if vOk == 0
    return
end

vSecondarySurface = cell(size(vPair2));

for f = 1:length(vPair2)
    vSecondarySurface{f} = vSurfacesList{vPair2(f)};
end

str=sprintf('Create Split Surfaces for %s', char(vPrimarySurface.GetName));
choice=questdlg(str, 'Create Split Surfaces',...
    'Analyze','Cancel','Analyze');
%Handle Response
switch choice
    case 'Analyze'
    case 'Cancel'
        return;
end

%%
%Get Image Data parameters
vImarisDataSet = vImarisApplication.GetDataSet.Clone;
vNumberOfChannels = vImarisDataSet.GetSizeC;
%vImarisDataSet.SetSizeC(vNumberOfChannels + 3);

vDataMin = [vImarisDataSet.GetExtendMinX, vImarisDataSet.GetExtendMinY, vImarisDataSet.GetExtendMinZ];
vDataMax = [vImarisDataSet.GetExtendMaxX, vImarisDataSet.GetExtendMaxY, vImarisDataSet.GetExtendMaxZ];
vDataSize = [vImarisDataSet.GetSizeX, vImarisDataSet.GetSizeY, vImarisDataSet.GetSizeZ];
aSizeT = vImarisApplication.GetDataSet.GetSizeT;
Xvoxelspacing = (vDataMax(1)-vDataMin(1))/vDataSize(1);
Yvoxelspacing = (vDataMax(2)-vDataMin(2))/vDataSize(2);
Zvoxelspacing = (vDataMax(3)-vDataMin(3))/vDataSize(3);


%Convert to 32bit
vFloatType = vImarisDataSet.GetType.eTypeFloat;
vImarisDataSet.SetType(vFloatType);



% Create a new channel where the result will be sent
vImarisDataSet.SetSizeC(vNumberOfChannels + length(vPair2));
%vImarisApplication.SetVisible(~vImarisApplication.GetVisible);

vMaskDataSetPrimary = vPrimarySurface.GetMask( ...
        vDataMin(1), vDataMin(2), vDataMin(3), ...
        vDataMax(1), vDataMax(2), vDataMax(3), ...
        vDataSize(1), vDataSize(2), vDataSize(3), 0);
    

vMaskDataSetSecondary = cell(size(vPair2));
for j = 1:length(vPair2)
    vMaskDataSetSecondary{j}= vSecondarySurface{j}.GetMask( ...
                  vDataMin(1), vDataMin(2), vDataMin(3), ...
                  vDataMax(1), vDataMax(2), vDataMax(3), ...
                  vDataSize(1), vDataSize(2), vDataSize(3), 0);
end

for vIndexZ = 1:vDataSize(3)
    %Mask for Primary surface (inside=1)
    vSlice = vMaskDataSetPrimary.GetDataSubVolumeAs1DArrayBytes(0,0,vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
    ch1 = double(vSlice);
       
    for k = 1:length(vPair2)
        %Mask for secondary Surface (inside=1)
        vSlice2=vMaskDataSetSecondary{k}.GetDataSubVolumeAs1DArrayBytes(0,0,vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
%        vSlice2 = vSlice2 == 1;% for inside surface only
        ch2 = double(vSlice2);
        vData=(ch1.*ch2);
        vImarisDataSet.SetDataSubVolumeAs1DArrayFloats(vData, ...
            0,0,vIndexZ-1,vNumberOfChannels+k-1,0,vDataSize(1),vDataSize(2),1);
    end
    
end

vImarisApplication.SetDataSet(vImarisDataSet);

%Create a new folder object for new surfaces
Newsurfaces = vImarisApplication.GetFactory;
result = Newsurfaces.CreateDataContainer;
result.SetName(sprintf('New Split Surfaces %s',vPrimarySurface.GetName));

for p = 1:length(vPair2)
    ip = vImarisApplication.GetImageProcessing;
    vNewSurface = ip.DetectSurfacesWithUpperThreshold(vImarisDataSet,[],vNumberOfChannels+p-1,0,0,true,false,0.99,true,false,1.01,'');
    vNewSurface.SetName(sprintf('%s-%s',vPrimarySurface.GetName,string(p)));
    vRGBA=[255,255,0, 0];%for yellow
    vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]); % combine different components (four bytes) into one integer
    vNewSurface.SetColorRGBA(vRGBA);
    result.AddChild(vNewSurface,-1);
    vImarisApplication.GetSurpassScene.AddChild(result,-1);
end

vImarisDataSet.Crop(0,vDataSize(1),0,vDataSize(2),0,vDataSize(3),0,vNumberOfChannels,0,aSizeT);

msgbox('New Split Surfaces Created!');
end
