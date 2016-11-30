function varargout = Sign_LanguageTranslator(varargin)
% SIGN_LANGUAGETRANSLATOR MATLAB code for Sign_LanguageTranslator.fig
%      SIGN_LANGUAGETRANSLATOR, by itself, creates a new SIGN_LANGUAGETRANSLATOR or raises the existing
%      singleton*.
%
%      H = SIGN_LANGUAGETRANSLATOR returns the handle to a new SIGN_LANGUAGETRANSLATOR or the handle to
%      the existing singleton*.
%
%      SIGN_LANGUAGETRANSLATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIGN_LANGUAGETRANSLATOR.M with the given input arguments.
%
%      SIGN_LANGUAGETRANSLATOR('Property','Value',...) creates a new SIGN_LANGUAGETRANSLATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Sign_LanguageTranslator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Sign_LanguageTranslator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Sign_LanguageTranslator

% Last Modified by GUIDE v2.5 26-Nov-2016 15:57:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Sign_LanguageTranslator_OpeningFcn, ...
                   'gui_OutputFcn',  @Sign_LanguageTranslator_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Sign_LanguageTranslator is made visible.
function Sign_LanguageTranslator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Sign_LanguageTranslator (see VARARGIN)

% Choose default command line output for Sign_LanguageTranslator
handles.output = hObject;
handles.datavalid = 0;
handles.data2valid = 0;
handles.extracted =0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Sign_LanguageTranslator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Sign_LanguageTranslator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in uploadsign.
function uploadsign_Callback(hObject, eventdata, handles)
% hObject    handle to uploadsign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Open a dialog box
[handles.file,handles.path] = uigetfile({'*.jpg';'*.tif';'*.*'}, 'Select a file from the directory', '../images/');

% Read the input image
handles.img = imread([handles.path, handles.file]);

% Modify the image
axes(handles.axes2);
imshow(handles.img);

handles.key = 1;

% Pre-Processing
[handles.imageop handles.data2valid] = preprocessing( handles.path, handles.file, handles.key);

% Complete
guidata(hObject, handles);

% --- Executes on button press in classify.
function classify_Callback(hObject, eventdata, handles)
% hObject    handle to classify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%Classify
handles=guidata(hObject);

if handles.data2valid~=0
    [ Orientation, Eccentricity, Width, Length, Fingers, Knuckles] = descriptor_calc( handles.imageop );
    handles.thisx = [ Orientation, Eccentricity, Width, Length, Fingers, Knuckles];
    
    % Compute the number of descriptor features.
    % NOTE: Overwrite existing for debugging
    numVectors = max(size(handles.thisx)) / 2;
    numVectors = 2;
    
    % Update the GUI
    set(handles.textOrientation,'string',Orientation);
    set(handles.textEccentricity,'string',Eccentricity);
    %set(handles.textFingers,'string',Fingers);
    set(handles.textWidth,'string',Width);
    set(handles.textLength,'string',Length);
    
    % Attempt to classify the input sign
    handles.label = classify( numVectors, handles.trainingFeatureVector, handles.C, handles.thisx ); % , handles.to, handles.t1);
    
    % Update the Output Image
    handles.img = getImage( '../imagesPristine/', handles.label );
    axes(handles.axes7);
    imshow(handles.img);
    
else
    warndlg('Error:No files uploaded')
end

% --- Executes on button press in loadimage.
function loadimage_Callback(hObject, eventdata, handles)
% hObject    handle to loadimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Overwriding user selection
%[handles.filename,handles.pathname] = uigetfile({'*.jpg';'*.tif';'*.*'}, 'Select a file from the directory');
%directory = handles.pathname;

directory = '../images/';
handles.filename = 0;
handles.key = 0;

[handles.imagestack handles.datavalid] = preprocessing( directory,handles.filename,handles.key);
% [ trainingFeatureVector, imageOrientation, imageEccentricity, imageWidth, imageLength, imageFingers, imageKnuckles] = descriptor_calc(handles. imagestack );
% x= [imageOrientation, imageEccentricity, imageWidth, imageLength, imageFingers, imageKnuckles];
%train x
%classify x

% if handles.filename == 0
%     handles.datavalid =1;
% else
%     handles.datavalid =0;
%     handles.img = imread([handles.pathname,handles.filename]);
% end
% axes(handles.axes1);
% imshow(handles.img);


 
% [File, Folder] = uigetfile('*.*', 'MultiSelect', 'on');
% 
% handles.img = cell(1, length(File));
% handles.axes = cell(1, length(File));
% for iFile = 1:length(File)
%     filename = fullfile(Folder, File{iFile});
%     image = imread(filename);
% %      axes(handles.axes(iFile));  
% %      imshow(image);
%     handles.img{iFile} = image;
% end
%  axes(handles.axes1);
%  imshow(handles.img{1});  
% main
  guidata(hObject, handles);

% --- Executes on button press in train.
function train_Callback(hObject, eventdata, handles)
% hObject    handle to train (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%train handles.x is feature vector

if handles.datavalid~=0 && handles.extracted ~=0
    handles = guidata(hObject);
    [to t1] = training(handles.x);
    handles.to = to;
    handles.t1 = t1;
elseif handles.datavalid == 0
    warndlg('Error:No files uploaded')
elseif handles.extracted ==0
    warndlg('Error:Please Extract Features first')
end
guidata(hObject, handles);


% --- Executes on button press in extractfeatures.
function extractfeatures_Callback(hObject, eventdata, handles)
% hObject    handle to extractfeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

if handles.datavalid~=0
    % Compute the training feature vectors
    [ imageOrientation, imageEccentricity, imageWidth, imageLength, imageFingers, imageKnuckles] = descriptor_calc( handles.imagestack );
    handles.x = [imageOrientation, imageEccentricity, imageWidth, imageLength, imageFingers, imageKnuckles];
    
    % Override the NUmber of feature vectors
    numVectors = max(size(handles.x)) / 2;
    numVectors = 2;

    [ handles.C ] = computeKMeansClusters( numVectors, handles.x );
    [ handles.trainingFeatureVector ] = computeTrainingVector(  numVectors, handles.C, handles.x );
    
    handles.extracted = 1;
    
    %plot ecentricity ish?
else
   warndlg('Error:No file uploaded') 
end

guidata(hObject, handles);
