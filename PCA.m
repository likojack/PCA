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

% Last Modified by GUIDE v2.5 03-Feb-2016 16:02:54

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

function plot_coord(axes3)
xlim = get(axes3, 'xlim');
ylim = get(axes3, 'ylim');
plot(axes3, xlim, [0 0], '--k');
plot(axes3, [0 0], ylim, '--k');


function [u,mad] = plot_pca(Data, axes4, popup)

xyz=ID2XYZ(deg2rad(Data.Inc),deg2rad(Data.Dec));
xyz=bsxfun(@times,xyz,Data.y);
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
[u,mad] = pmag_pca(x,y,z);
[coeff,score]=princomp(xyz(:,:));

hat=[min(score(:,1));max(score(:,1))]*coeff(:,1)';
%hat=bsxfun(@plus,hat,mean(xyz(id1,:)));
hat=bsxfun(@plus,hat,mean(xyz(:,:)));

contents = get(popup,'String');
content = contents{get(popup,'Value')};
switch content
    case 'V vs N'
        plot(axes4, xyz(:,2),xyz(:,1),'ok','markerfacecolor','k');
        hold(axes4,'on');
        plot(axes4, xyz(:,1),-xyz(:,3),'sk','markerfacecolor','w');
        plot(axes4, hat(:,2),hat(:,1),'r','linewidth',1.5);
        plot(axes4, hat(:,1),-hat(:,3),'r','linewidth',1.5);
        set(axes4,'Dataaspectratio', [1 1 1]);
        plot_coord(axes4);
        hold(axes4, 'off');
    case 'V vs E'
        plot(axes4, xyz(:,2),xyz(:,1),'ok','markerfacecolor','k');
        hold(axes4,'on');
        plot(axes4, xyz(:,2),-xyz(:,3),'sk','markerfacecolor','w');
        plot(axes4, hat(:,2),hat(:,1),'r','linewidth',1.5);
        plot(axes4, hat(:,2),-hat(:,3),'r','linewidth',1.5);
        set(axes4,'Dataaspectratio', [1 1 1]);
        plot_coord(axes4);
        hold(axes4,'off');
     case 'V vs H'
        plot(axes4, xyz(:,2),xyz(:,1),'ok','markerfacecolor','k');
        hold(axes4,'on');
        plot(axes4, sqrt(xyz(:,1).^2 + xyz(:,2).^2),-xyz(:,3),'sk','markerfacecolor','w');
        plot(axes4, hat(:,2),hat(:,1),'r','linewidth',1.5);
        plot(axes4, sqrt(hat(:,1).^2 + hat(:,2).^2),-hat(:,3),'r','linewidth',1.5);
        set(axes4,'Dataaspectratio', [1 1 1]);
        plot_coord(axes4);
        hold(axes4,'off');
end


function plot_third(dat, Dat, axes3, popup_direction)
color_xyz = ID2XYZ(deg2rad(Dat.Inc),deg2rad(Dat.Dec));
color_xyz = bsxfun(@times, color_xyz, Dat.y);

xyz = ID2XYZ(deg2rad(dat.Inc),deg2rad(dat.Dec));
xyz = bsxfun(@times, xyz, dat.y);

contents = get(popup_direction,'String');
content = contents{get(popup_direction,'Value')};
switch content
    case 'V vs N'
        plot(axes3, xyz(:,2),xyz(:,1),'ok-','markerfacecolor','k');
        hold(axes3,'on');
        plot(axes3, xyz(:,1),-xyz(:,3),'sk-','markerfacecolor','w')
        plot(axes3, color_xyz(:,2),color_xyz(:,1),'*', 'color', 'r');
        plot(axes3, color_xyz(:,1),-color_xyz(:,3),'*', 'color', 'b');
        set(axes3,'Dataaspectratio', [1 1 1]);
        plot_coord(axes3);
        hold(axes3,'off');
    case 'V vs E'
        plot(axes3, xyz(:,2),xyz(:,1),'ok-','markerfacecolor','k');
        hold(axes3,'on');
        plot(axes3, xyz(:,2),-xyz(:,3),'sk-','markerfacecolor','w');
        plot(axes3, color_xyz(:,2),color_xyz(:,1),'*', 'color', 'r');
        plot(axes3, color_xyz(:,2),-color_xyz(:,3),'*', 'color', 'b');
        plot_coord(axes3);
        set(axes3,'Dataaspectratio', [1 1 1]);
        hold(axes3,'off');   
     case 'V vs H'
        plot(axes3, xyz(:,2),xyz(:,1),'ok-','markerfacecolor','k');
        hold(axes3,'on');
        plot(axes3, sqrt(xyz(:,1).^2 + xyz(:,2).^2),-xyz(:,3),'sk-','markerfacecolor','w')
        plot(axes3, color_xyz(:,2),color_xyz(:,1),'*', 'color', 'r');
        plot(axes3, sqrt(color_xyz(:,1).^2 + color_xyz(:,2).^2),-color_xyz(:,3),'*', 'color', 'b');
        set(axes3,'Dataaspectratio', [1 1 1]);
        plot_coord(axes3);
        hold(axes3,'off');    
end

function plot_first(x,y,axes1)
plot(axes1, x, y, '*-');


function plot_second(dat, axes2)
bx = findobj('Tag', 'axes2');
cla(axes2,'reset');
hold(axes2,'on');
I0=[0:15:90]';
for i=1:numel(I0)
    D1=[0:1:360]';
    I1=ones(size(D1)).*I0(i);
    XYZtick=ID2XYZ(deg2rad(I1),deg2rad(D1));
    Stick=pmag_splot(XYZtick);
    plot(axes2, Stick(:,1),Stick(:,2),'k');
end
set(axes2,'box','off','visible','off');
XYZ = ID2XYZ(dat.Inc, dat.Dec);
XYZ = bsxfun(@times, dat.y,XYZ);
% XYZ = [dat.x_pos, dat.y_pos, dat.z_pos];

XYZ=bsxfun(@rdivide,XYZ,sqrt(sum(XYZ.^2,2)));
Sxy = pmag_splot(XYZ);
plot(axes2, Sxy(:,1),Sxy(:,2),'or');
hold(axes2,'off');


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

%%get the seleected data index
index_selected = get(hObject, 'Value');
%%get the selected data
h = findobj('Tag', 'Load_dataset');
dat = get(h,'UserData');
x = dat.x([index_selected]);
y = dat.y([index_selected]);
Inc = dat.Inc([index_selected]);
Dec = dat.Dec([index_selected]);
Data = struct('x', x, 'y', y, 'Inc', Inc, 'Dec', Dec);
setappdata(hObject, 'data', Data);
%%plot first diagram
plot_first(x,y,handles.axes1);
%%plot second diagram
plot_second(Data, handles.axes2);


%%plot the third diagram

plot_third(dat, Data, handles.axes3, handles.popup_direction);


%%project data from different directions
[u,mad] = plot_pca(Data, handles.axes4, handles.popup_direction);

set(handles.MAD, 'String', mad);
[I,D] = XYZ2ID(u');
set(handles.Inclination, 'String' ,rad2deg(I));
set(handles.Declination, 'String', rad2deg(D));


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

% xin=[0.00 9.2830e-05 339.9 57.9
% 2.50 7.5820e-05 325.7 49.1
% 5.00 6.2920e-05 321.3 45.9
% 10.00 5.2090e-05 314.8 41.7
% 15.00 4.4550e-05 310.3 38.7
% 20.00 3.9540e-05 305.0 37.0
% 30.00 3.2570e-05 303.9 34.7
% 40.00 2.5670e-05 303.0 32.3
% 50.00 2.2520e-05 303.6 32.4
% 60.00 1.9820e-05 299.8 30.8
% 70.00 1.3890e-05 292.5 31.0
% 80.00 1.2570e-05 297.0 25.6
% 90.00 0.5030e-05 299.3 11.3];
% x = xin(:,1);
% y = xin(:,2);
% Dec = xin(:,3);
% Inc = xin(:,4);

%x,y here show the values of first diagram.
x = rawdata1.AFD0x28mT0x29;
y = rawdata1.Intensity0x28tray0x29;

%x_pos, y_pos, z_pos here show the position in 3d space.
x_pos = rawdata1.XtrayC;
y_pos = rawdata1.YtrayC;
z_pos = rawdata1.ZtrayC;

Dec = rawdata1.Declination0x28tray0x29;
Inc = rawdata1.Inclination0x28tray0x29;

% Dat = struct('name', name, 'x', x, 'y', y, 'x_pos', x_pos, 'y_pos', y_pos, 'z_pos', z_pos, 'Inc', Inc, 'Dec', Dec);
Dat = struct('name', name, 'x', x, 'y', y,'Inc', Inc, 'Dec', Dec);

set(hObject, 'UserData', Dat);
n = size(Dec);
set(handles.listbox1,'String', x);

% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listbox1, 'String', '');
cla(handles.axes1,'reset');
cla(handles.axes2,'reset');
cla(handles.axes3,'reset');
cla(handles.axes4,'reset');

function Inclination_Callback(hObject, eventdata, handles)
% hObject    handle to Inclination (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of Inclination as text
%        str2double(get(hObject,'String')) returns contents of Inclination as a double


% --- Executes during object creation, after setting all properties.
function Inclination_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Inclination (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Declination_Callback(hObject, eventdata, handles)
% hObject    handle to Declination (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of Declination as text
%        str2double(get(hObject,'String')) returns contents of Declination as a double


% --- Executes during object creation, after setting all properties.
function Declination_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Declination (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MAD_Callback(hObject, eventdata, handles)
% hObject    handle to MAD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MAD as text
%        str2double(get(hObject,'String')) returns contents of MAD as a double


% --- Executes during object creation, after setting all properties.
function MAD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MAD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_direction.
function popup_direction_Callback(hObject, eventdata, handles)
% hObject    handle to popup_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_direction contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_direction

%%decide which data source it is coming from

data_list = getappdata(handles.listbox1,'data');
plot_first(data_list.x, data_list.y, handles.axes1);

plot_second(data_list, handles.axes2);
dat_h = findobj('Tag', 'Load_dataset');
data_run = get(dat_h,'UserData');    
plot_third(data_run, data_list, handles.axes3, handles.popup_direction);
plot_pca(data_list, handles.axes4, handles.popup_direction);

% --- Executes during object creation, after setting all properties.
function popup_direction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MAP_I_Callback(hObject, eventdata, handles)
% hObject    handle to MAP_I (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MAP_I as text
%        str2double(get(hObject,'String')) returns contents of MAP_I as a double


% --- Executes during object creation, after setting all properties.
function MAP_I_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MAP_I (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mD0_Callback(hObject, eventdata, handles)
% hObject    handle to mD0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mD0 as text
%        str2double(get(hObject,'String')) returns contents of mD0 as a double


% --- Executes during object creation, after setting all properties.
function mD0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mD0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ip1_Callback(hObject, eventdata, handles)
% hObject    handle to Ip1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ip1 as text
%        str2double(get(hObject,'String')) returns contents of Ip1 as a double


% --- Executes during object creation, after setting all properties.
function Ip1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ip1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MAP_D_Callback(hObject, eventdata, handles)
% hObject    handle to MAP_D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MAP_D as text
%        str2double(get(hObject,'String')) returns contents of MAP_D as a double


% --- Executes during object creation, after setting all properties.
function MAP_D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MAP_D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theta_Callback(hObject, eventdata, handles)
% hObject    handle to theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theta as text
%        str2double(get(hObject,'String')) returns contents of theta as a double


% --- Executes during object creation, after setting all properties.
function theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Dp_Callback(hObject, eventdata, handles)
% hObject    handle to Dp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dp as text
%        str2double(get(hObject,'String')) returns contents of Dp as a double


% --- Executes during object creation, after setting all properties.
function Dp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Bayesian.
function Bayesian_Callback(hObject, eventdata, handles)
% hObject    handle to Bayesian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get data point structure
data = getappdata(handles.listbox1,'data');
xyz=ID2XYZ(deg2rad(data.Inc),deg2rad(data.Dec));
xyz=bsxfun(@times,xyz,data.y);

p = str2double(get(handles.confidence,'String'));
data_list = getappdata(handles.listbox1,'data');
inc = deg2rad(data_list.Inc);
dec = deg2rad(data_list.Dec);
XYZ = ID2XYZ(inc,dec);
XYZ = bsxfun(@times, data_list.y,XYZ);

%show up wait window and disable user input
h = msgbox('This could take a while. Patience');
set(handles.listbox1,'Enable','off');
set(handles.Load_dataset,'Enable','off');
set(handles.clear,'Enable','off');
set(handles.Bayesian,'Enable','off');
set(handles.popup_direction,'Enable','off');
set(hObject,'Interruptible','off');



[MAP_I,MAP_D,theta95,MD0,Ip,Dp,VN,VE,VH,Emu] = pmag_bpca(XYZ(:,1),XYZ(:,2),XYZ(:,3),p);
MAP_I= deg2rad(MAP_I);
MAP_D = deg2rad(MAP_D);
map_xyz = ID2XYZ(MAP_I, MAP_D);

min_pt = min(VN.H{1}(:,2))*map_xyz;
max_pt = max(VN.H{1}(:,2))*map_xyz;
hat=bsxfun(@plus,[min_pt;max_pt],Emu');

%scale the VN.H line
k1 = (hat(1,1)-hat(2,1))/(hat(1,2)-hat(2,2));
b1 = hat(1,1)-k1*hat(1,2);
x1_min = (min(VN.H{1}(:,2)) - b1)/k1;
x1_max = (max(VN.H{1}(:,2)) - b1)/k1;

%scale the VN.V line
k2 = (-hat(1,3)+hat(2,3))/(hat(1,1)-hat(2,1));
b2 = -hat(1,3)-k2*hat(1,1);
x2_min = (min(VN.V{1}(:,2)) - b2)/k2;
x2_max = (max(VN.V{1}(:,2)) - b2)/k2;

%scale VE.H line
k3 = k1;
b3 = b1;
x3_min = (min(VE.H{1}(:,2)) - b3)/k3;
x3_max = (max(VE.H{1}(:,2)) - b3)/k3;

%scale VE.V line
k4 = (-hat(1,3)+hat(2,3))/(hat(1,2)-hat(2,2));
b4 = -hat(1,3) - k4*hat(1,2);
x4_min = (min(VE.V{1}(:,2)) - b4)/k4;
x4_max = (max(VE.V{1}(:,2)) - b4)/k4;

%scale VH.H line
k5 = k1;
b5 = b1;
x5_min = (min(VH.H{1}(:,2)) - b5)/k5;
x5_max = (max(VH.H{1}(:,2)) - b5)/k5;

%scale VH.V line
k6 = (-hat(1,3)+hat(2,3))/(sqrt(hat(1,1).^2+hat(1,2).^2) - sqrt(hat(2,1).^2+hat(2,2).^2));
b6 = -hat(1,3)-k6*sqrt(hat(1,1).^2+hat(1,2).^2);
x6_min = (min(VH.V{1}(:,2)) - b6)/k6;
x6_max = (max(VH.V{1}(:,2)) - b6)/k6;

close(h);
set(handles.listbox1,'Enable','on');
set(handles.Load_dataset,'Enable','on');
set(handles.clear,'Enable','on');
set(handles.Bayesian,'Enable','on');
set(handles.popup_direction,'Enable','on');

contents = get(handles.popup_direction,'String');
content = contents{get(handles.popup_direction,'Value')};
switch content
    case 'V vs N'
        plot(handles.axes4, VN.H{1}(:,1), VN.H{1}(:,2),'--k');
        hold(handles.axes4,'on');
        plot(handles.axes4,[x1_min;x1_max],[min(VN.H{1}(:,2));max(VN.H{1}(:,2))]);
        plot(handles.axes4,VN.H{2}(:,1), VN.H{2}(:,2),'--k');
        plot(handles.axes4,xyz(:,2),xyz(:,1),'ok','markerfacecolor','k');

        plot(handles.axes4,VN.V{1}(:,1), VN.V{1}(:,2),'--r');
        plot(handles.axes4,VN.V{2}(:,1), VN.V{2}(:,2),'--r');
        plot(handles.axes4,[x2_min;x2_max],[min(VN.V{1}(:,2));max(VN.V{1}(:,2))]);
        plot(handles.axes4,xyz(:,1),-xyz(:,3),'sk','markerfacecolor','w');
        set(handles.axes4,'Dataaspectratio', [1 1 1]);
        plot_coord(handles.axes4);
        hold(handles.axes2,'off');
    case 'V vs E'
        plot(handles.axes4, VE.H{1}(:,1), VE.H{1}(:,2),'--k');
        hold(handles.axes4,'on');
        plot(handles.axes4,VE.H{2}(:,1),VE.H{2}(:,2),'--k');
        plot(handles.axes4,[x3_min;x3_max],[min(VE.H{1}(:,2));max(VE.H{1}(:,2))]);
        plot(handles.axes4,xyz(:,2),xyz(:,1),'ok','markerfacecolor','k');

        plot(handles.axes4,VE.V{1}(:,1), VE.V{1}(:,2),'--r');
        plot(handles.axes4,VE.V{2}(:,1), VE.V{2}(:,2),'--r');
        plot(handles.axes4,[x4_min;x4_max],[min(VE.V{1}(:,2));max(VE.V{1}(:,2))]);
        plot(handles.axes4,xyz(:,2),-xyz(:,3),'sk','markerfacecolor','w');
        set(handles.axes4,'Dataaspectratio', [1 1 1]);
        plot_coord(handles.axes4);
        hold(handles.axes4,'off');        
    case 'V vs H'
        plot(handles.axes4, VH.H{1}(:,1), VH.H{1}(:,2),'--k');
        hold(handles.axes4,'on');
        plot(handles.axes4,VH.H{2}(:,1), VH.H{2}(:,2),'--k');
        plot(handles.axes4,[x5_min;x5_max],[min(VN.H{1}(:,2));max(VN.H{1}(:,2))]);
        plot(handles.axes4,xyz(:,2),xyz(:,1),'ok','markerfacecolor','k');
        % 
        plot(handles.axes4,VH.V{1}(:,1), VH.V{1}(:,2),'--r');
        plot(handles.axes4,VH.V{2}(:,1), VH.V{2}(:,2),'--r');
        plot(handles.axes4,[x6_min;x6_max],[min(VH.V{1}(:,2));max(VH.V{1}(:,2))]);
        plot(handles.axes4, sqrt(xyz(:,1).^2 + xyz(:,2).^2),-xyz(:,3),'sk','markerfacecolor','w');
        set(handles.axes4,'Dataaspectratio', [1 1 1]);
        plot_coord(handles.axes4);
        hold(handles.axes4,'off');

end



set(handles.MAP_I,'String',MAP_I);
set(handles.MAP_D,'String',MAP_D);
set(handles.Ip1,'String',Ip(1));
set(handles.Dp1,'String',Dp(1));
set(handles.Ip2,'String',Ip(2));
set(handles.Dp2,'String',Dp(2));
set(handles.mD0,'String',MD0);
set(handles.theta,'String',theta95);



function Ip2_Callback(hObject, eventdata, handles)
% hObject    handle to Ip2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ip2 as text
%        str2double(get(hObject,'String')) returns contents of Ip2 as a double


% --- Executes during object creation, after setting all properties.
function Ip2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ip2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Dp1_Callback(hObject, eventdata, handles)
% hObject    handle to Dp1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dp1 as text
%        str2double(get(hObject,'String')) returns contents of Dp1 as a double


% --- Executes during object creation, after setting all properties.
function Dp1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dp1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Dp2_Callback(hObject, eventdata, handles)
% hObject    handle to Dp2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dp2 as text
%        str2double(get(hObject,'String')) returns contents of Dp2 as a double


% --- Executes during object creation, after setting all properties.
function Dp2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dp2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function confidence_Callback(hObject, eventdata, handles)
% hObject    handle to confidence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of confidence as text
%        str2double(get(hObject,'String')) returns contents of confidence as a double


% --- Executes during object creation, after setting all properties.
function confidence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to confidence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
