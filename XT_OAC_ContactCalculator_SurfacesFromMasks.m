%Create Surfaces from Masks

%Written by Olivia Creasey, PhD Student.
%This version written September 14-15, 2023. 
%Based off of previous "Surfaces from Multicut" XTensions written
%August 2019, Edit May 2023 and July 2023 and September 2023

%
%<CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Create Surfaces from Masks" icon="Matlab">
%          <Command>MatlabXT::XT_OAC_ContactCalculator_SurfacesFromMasks(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Create Surfaces from Masks" icon="Matlab">
%            <Command>MatlabXT::XT_OAC_ContactCalculator_SurfacesFromMasks(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%Description
%This XTension creates individual surface objects from a specified
%image channel. The image channel in this case is a channel containing 
%masks, specifically the output of an ilastik Multicut workflow, in which 
%all pixels in the image volume corresponding to individual objects (cells
%etc) are annotated as having the same intensity value. Each object has its
%own intensity value. Masks for an indefinite number of objects can be
%contained in the same image channel. This XTension was written as part of 
%the suite of XTensions written for the Contact Calculator analysis but can
%be used outside of the Contact Calculator analysis.
%Written for use with Imaris 9.5.0

function XT_OAC_ContactCalculator_SurfacesFromMasks(aImarisApplicationID)

%connect to Imaris interface
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

%Get Image Data parameters
vImarisDataSet = vImarisApplication.GetDataSet.Clone;
vNumberOfChannels = vImarisDataSet.GetSizeC;
vDataSize = [vImarisDataSet.GetSizeX, vImarisDataSet.GetSizeY, vImarisDataSet.GetSizeZ];

vChannelsList = 1:vNumberOfChannels;

%Create Dialog box and allow user to choose the Channel containing the
%Multicut output
[vSegmentationChannel, vOk] = listdlg('ListString',string(vChannelsList),'SelectionMode','single',...
    'ListSize',[250 150],'Name','Select Multicut Channel','InitialValue',1, ...
    'PromptString',{'Please select the Channel containing the Multicut Output'});

if vOk == 0
    return
end

str=sprintf('Create Surfaces from Channel %s?\n', string(vChannelsList(vSegmentationChannel)));
choice=questdlg(str, 'Create surfaces',...
    'Analyze','Cancel','Analyze');
%Handle Response
switch choice
    case 'Analyze'
    case 'Cancel'
        return;
end


%Only some intensity values are used, so we need to know what they are
    %Early versions of the XTension attempted to use the Imaris internal
    %IDataSet.GetChannelRangeMin and IDataSet.GetChannelRangeMax functions 
    %and it didn't work for unknown reasons.
vChannelIntensityList = [];
%Java can't handle large data sets, so I have to do this by slice
for vIndexZ = 1:vDataSize(3)
    vSliceIntensityList = vImarisDataSet.GetDataSubVolumeAs1DArrayFloats(0,0, vIndexZ-1, vSegmentationChannel-1 ,0, vDataSize(1), vDataSize(2), 1);
    vChannelIntensityList = unique([reshape(vChannelIntensityList,1,[]), reshape(vSliceIntensityList,1,[])]);  
end

vChannelIntensityCount = length(vChannelIntensityList);

%Create a new folder object for new surfaces
Newsurfaces = vImarisApplication.GetFactory;
result = Newsurfaces.CreateDataContainer;
result.SetName(sprintf('Surfaces from Multicut Channel %s',string(vSegmentationChannel)));

%Create a waitbar so that the user can monitor progress even with the
%window made invisible
vProgressDisplay = waitbar(0, 'Producing Surface Objects');
vImarisApplication.SetVisible(0); %This speeds the process up some AND 
%reduces the amount of processing power and RAM needed


%Create new surfaces, producing only those surfaces corresponding to
%intensity values actually in the Multicut channel
ip = vImarisApplication.GetImageProcessing;
vRGBA=[255,255,255, 0];%for yellow
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]); % combine different components (four bytes) into one integer
for b = 1:vChannelIntensityCount
    
    UpperThreshold = vChannelIntensityList(b)+0.5;
    LowerThreshold = vChannelIntensityList(b)-0.5;
    vNewSurface = ip.DetectSurfacesWithUpperThreshold(vImarisDataSet,[],vSegmentationChannel-1,0,0,true,false,LowerThreshold,true,false,UpperThreshold,'');
    vNewSurface.SetName(sprintf('Surface %s',string(vChannelIntensityList(b))));
    vNewSurface.SetColorRGBA(vRGBA);
    %Add new surface to Surpass Scene
    result.AddChild(vNewSurface, -1);
    vImarisApplication.GetSurpassScene.AddChild(result, -1);

    waitbar(b/vChannelIntensityCount,vProgressDisplay);
    
end
    
close(vProgressDisplay);
vImarisApplication.SetVisible(1);
msgbox('Multicut Surfaces Generated!');

end
