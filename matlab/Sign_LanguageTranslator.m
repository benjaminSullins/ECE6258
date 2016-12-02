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
handles.trainingdata=0;
currentFolder = cd;

%fileExisting  = (exist(fullfile(currentFolder, 'File'), 'file') == 2);
file ='trainingFeatureVectorDefault.mat';
% file2 ='CDefault.mat';
if exist(fullfile(currentFolder, file), 'file') == 2
    disp('file found')
    aa= load('trainingFeatureVectorDefault.mat', 'CDefault');
    bb= load('trainingFeatureVectorDefault.mat', 'trainingFeatureVectorDefault');
    handles.trainingFeatureVector = bb.trainingFeatureVectorDefault
    disp('file2 found')
    disp(handles.trainingFeatureVector)
    C = aa.CDefault;
    
     handles.C = C;
     disp(handles.trainingFeatureVector)
% %     handles.C = num2cell(cluster,[1,6]);
%      handles.C = cluster;
    handles.trainingdata=1;
    
end

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

if handles.data2valid~=0 &&(handles.extracted == 1 ||handles.trainingdata ==1)
   
            disp('testing')
            if(handles.extracted ==0)
            if isempty(handles.trainingFeatureVector)|| isempty(handles.C)
                errordlg('Empty Set:Extract Feature Vectors First');
                return
            else
                warndlg('Using Default Training Vectors');
               
            end
            end
       
    
    %[ Orientation, Eccentricity, Width, Length, Fingers, Knuckles] = descriptor_calc( handles.imageop );
    %handles.thisx = [ Orientation, Eccentricity, Width, Length, Fingers, Knuckles];
 [ Orientation, Eccentricity, Width, Length, Fingers, Knuckles,Fourier_mean,Fourier_max,Fourier_sigma,Fourier_min,Fourier_dc,Fourier_first] = descriptor_calc( handles.imageop);
    %handles.thisx =[ Orientation, Eccentricity, Width, Length, Fingers, Knuckles,real(Fourier_mean),imag(Fourier_mean),real(Fourier_max),imag(Fourier_max),real(Fourier_sigma),imag(Fourier_sigma),real(Fourier_min),imag(Fourier_min),real(Fourier_dc),imag(Fourier_dc),real(Fourier_first),imag(Fourier_first), Fingers, Knuckles];

    handles.thisx =[ Orientation, Eccentricity, Width, Length,real(Fourier_max),imag(Fourier_max),real(Fourier_min),imag(Fourier_min),real(Fourier_dc),imag(Fourier_dc),real(Fourier_first),imag(Fourier_first), Fingers, Knuckles];



    % Compute the number of descriptor features.
    % NOTE: Overwrite existing for debugging
    numVectors = max(size(handles.thisx)) / 2;
    %numVectors = 2;
    %numVectors = 8;

    numVectors = 6;

    
    % Update the GUI
    set(handles.textOrientation,'string',Orientation);
    set(handles.textEccentricity,'string',Eccentricity);
    %set(handles.textFingers,'string',Fingers);
    set(handles.textWidth,'string',Width);
    set(handles.textLength,'string',Length);
    
    % Attempt to classify the input sign
    [handles.label handles.corrPlot] = classify( numVectors, handles.trainingFeatureVector, handles.C, handles.thisx ); % , handles.to, handles.t1);
    
    % Update the Output Image
    handles.img = getImage( '../imagesPristine/', handles.label );

%     axes(handles.axes7);
%     imshow(handles.img);

    % Plot the correlation
    %handles.disp = plot(handles.corrPlot,'Parent',handles.axes4);
    axes(handles.axes3);
     plot( handles.corrPlot,'k*','MarkerSize',5);
     title('Correlation')

    
    % Clear the axes to prevent ghost images overlaying each other
    cla(handles.axes7);
    
    % Update axes with the pristine image
    axes(handles.axes7);
    imshow(handles.img);
%     handles.disp = imshow(handles.img,'Parent',handles.axes7);

    
elseif handles.data2valid==0  
    warndlg('Error:No files uploaded');
elseif handles.extracted == 0
     warndlg('Error:Please Extract features first');   
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

  guidata(hObject, handles);

% --- Executes on button press in train.
function train_Callback(hObject, eventdata, handles)
% hObject    handle to train (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%train handles.x is feature vector

if handles.datavalid~=0 && handles.extracted ~=0
%     handles = guidata(hObject);
%     [to t1] = training(handles.x);
%     handles.to = to;
%     handles.t1 = t1;
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
    [ imageOrientation, imageEccentricity, imageWidth, imageLength, imageFingers, imageKnuckles,Fourier_mean,Fourier_max,Fourier_sigma,Fourier_min,Fourier_dc,Fourier_first] = descriptor_calc( handles.imagestack );
    handles.x =[ imageOrientation, imageEccentricity, imageWidth, imageLength,real(Fourier_max),imag(Fourier_max),real(Fourier_min),imag(Fourier_min),real(Fourier_dc),imag(Fourier_dc),real(Fourier_first),imag(Fourier_first), imageFingers, imageKnuckles];


    
    % Override the NUmber of feature vectors
    numVectors = max(size(handles.x)) / 2;
    %numVectors = 2;
    %numVectors = 8;
    numVectors = 6;


    [ C , handles.XGrid,handles.idx2Region,handles.X] = computeKMeansClusters( numVectors, handles.x );
    handles.C = C;
    [ trainingFeatureVector ] = computeTrainingVector(  numVectors, handles.C, handles.x );
    
      
      handles.trainingFeatureVector = trainingFeatureVector;
%     save('CDefault');
%     save('trainingFeatureVector');

    % Assigns each node in the grid to the closest centroid    
    
    
    gscatter(handles.XGrid(:,1),handles.XGrid(:,2),handles.idx2Region,...
        [0,0.75,0.75;0.75,0,0.75;0.75,0.75,0],'..');
    hold on;
    axes(handles.axes1);
    plot(handles.X(:,1),handles.X(:,2),'k*','MarkerSize',5);
    title 'Width vs. Length';
    xlabel 'Width (pixels)';
    ylabel 'Length (pixels)';
    legend('Region 1','Region 2','Region 3','Region 4','Region 5','Data','Location','SouthEast');
    hold off;
    
    
    % Plot Frequency Domain Derived Descriptors
axes(handles.axes4);
plot((real(Fourier_max)),'*')
hold on
plot(real(Fourier_min),'o')
hold on
plot(real(Fourier_dc),'+')
hold on
plot(real(Fourier_first),'v')
plot((imag(Fourier_max)),'*')
hold on
plot(imag(Fourier_min),'o')
hold on
plot(imag(Fourier_dc),'+')
hold on
plot(imag(Fourier_first),'v');
legend('x-max','x-min','x-Zero-frequency','x-First Harmonic','y-max','y-min','y-Zero-frequency','y-First Harmonic')
title('Frequency Domain Descriptors');
    
    handles.extracted = 1;
    
    %plot ecentricity ish?
else
   warndlg('Error:No file uploaded') 
end

guidata(hObject, handles);
