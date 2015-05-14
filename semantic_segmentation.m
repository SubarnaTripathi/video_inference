function varargout = semantic_segmentation(varargin)
% SEMANTIC_SEGMENTATION MATLAB code for semantic_segmentation.fig
%      SEMANTIC_SEGMENTATION, by itself, creates a new SEMANTIC_SEGMENTATION or raises the existing
%      singleton*.
%
%      H = SEMANTIC_SEGMENTATION returns the handle to a new SEMANTIC_SEGMENTATION or the handle to
%      the existing singleton*.
%
%      SEMANTIC_SEGMENTATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEMANTIC_SEGMENTATION.M with the given input arguments.
%
%      SEMANTIC_SEGMENTATION('Property','Value',...) creates a new SEMANTIC_SEGMENTATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before semantic_segmentation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to semantic_segmentation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help semantic_segmentation

% Last Modified by GUIDE v2.5 22-Jul-2014 17:04:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @semantic_segmentation_OpeningFcn, ...
                   'gui_OutputFcn',  @semantic_segmentation_OutputFcn, ...
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


% --- Executes just before semantic_segmentation is made visible.
function semantic_segmentation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to semantic_segmentation (see VARARGIN)

% Choose default command line output for semantic_segmentation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);




% UIWAIT makes semantic_segmentation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = semantic_segmentation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in semantic_segmentation.
function semantic_segmentation_Callback(hObject, eventdata, handles)
% hObject    handle to semantic_segmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes3);
imshow('CamVid\legends.png');

%%% RUN AS A WHOLE
%system('x64\Debug\DenseHO.exe');

base_num = 63;
base_input_str = 'CamVid\Result\CImage\001TP_';
base_output_str = 'CamVid\Result\Crf\001TP_';

% base_input_str = 'CImage\001TP_';
% base_output_str = 'Crf\001TP_';


for im_idx = 0:5 %37  
    handles.frame_id = base_num+im_idx;
    
    str = sprintf('Frame # %d', handles.frame_id);
    set(findobj('Tag','frame_num'),'String',str);
    
    %%%%% call one test frame at a time
    in_cmd = sprintf('x64\\Debug\\DenseHO_demo.exe 001TP_%03d',base_num+im_idx);
    
    
    %%% run as individual file
    [status,cmdout]= system(in_cmd);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    in_img = sprintf('%s%03d.jpg', base_input_str, base_num+im_idx);
    out_img = sprintf('%s%03d.png', base_output_str, base_num+im_idx);

    % Update handles structure
    guidata(hObject, handles);

    axes(handles.axes1);
    %imshow('sample_image1.jpg');
    imshow(in_img);

    axes(handles.axes2);
    imshow(out_img);
    
         
%     GUI(set(handles.frame_num, 'String', da ));
%     guidata(hObject, handles);
    
    %pause(1.0);
end




% --- Executes on button press in stop_process.
function stop_process_Callback(hObject, eventdata, handles)
% hObject    handle to stop_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pause


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function frame_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% str = sprintf('frame num = %d', handles.frame_id);
% GUI(set(handles.frame_num, 'String', str));


% --- Executes during object deletion, before destroying properties.
function frame_num_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to frame_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
