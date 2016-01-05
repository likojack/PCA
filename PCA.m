function varargout = PCA(varargin)
% PCA MATLAB code for PCA.fig
%      PCA, by itself, creates a new PCA or raises the existing
%      singleton*.
%
%      H = PCA returns the handle to a new PCA or the handle to
%      the existing singleton*.
%
%      PCA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PCA.M with the given input arguments.
%
%      PCA('Property','Value',...) creates a new PCA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PCA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PCA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PCA

% Last Modified by GUIDE v2.5 05-Jan-2016 13:30:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PCA_OpeningFcn, ...
                   'gui_OutputFcn',  @PCA_OutputFcn, ...
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


% --- Executes just before PCA is made visible.
function PCA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PCA (see VARARGIN)

% Choose default command line output for PCA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PCA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PCA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
index_selected = get(hObject, 'Value');
h = findobj('Tag', 'Load_dataset');
dat = get(h,'UserData');
x = dat.x([index_selected]);
y = dat.y([index_selected]);
plot(handles.axes1, x, y, '*-');


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in Bayesian.
function Bayesian_Callback(hObject, eventdata, handles)
% hObject    handle to Bayesian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Load_dataset.
function Load_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to Load_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%load data here, what data format?
[filename1,filepath1]=uigetfile({'*.*','All Files'},...
  'Select Data File 1');
cd(filepath1);
rawdata1 = tdfread(filename1);
name = rawdata1.Sample;
%x,y here show the values of first diagram.
x = rawdata1.AFD0x28mT0x29;
y = rawdata1.Intensity0x28tray0x29;

%x_pos, y_pos, z_pos here show the position in 3d space.
x_pos = rawdata1.XtrayC;
y_pos = rawdata1.YtrayC;
z_pos = rawdata1.ZtrayC;

Dec = rawdata1.Declination0x28tray0x29;
Inc = rawdata1.Inclination0x28tray0x29;
Dat = struct('name', name, 'x', x, 'y', y, 'Dec', Dec, 'Inc', Inc);
set(hObject, 'UserData', Dat);
n = size(Dec);
set(handles.listbox1,'String', x);


% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ax = findobj('Tag', 'axes1');
dat_h = findobj('Tag', 'Load_dataset');
dat = get(dat_h,'UserData');
plot(ax,dat.x,dat.y, '*-');



% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listbox1, 'String', '');
cla(handles.axes1,'reset');
cla(handles.axes2,'reset');
cla(handles.axes3,'reset');
