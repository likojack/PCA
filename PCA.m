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

% Last Modified by GUIDE v2.5 11-Jan-2016 09:33:06

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

function plot_pca(xyz,hat, axes3, axes4, popup,id1)
contents = get(popup,'String');
content = contents{get(popup,'Value')}
switch content
    case 'V vs N'
        plot(axes3, xyz(:,2),xyz(:,1),'ok-','markerfacecolor','k');
        hold(axes3,'on');
        plot(axes3, xyz(:,1),-xyz(:,3),'sk-','markerfacecolor','w')
        hold(axes3,'off');

        plot(axes4, xyz(id1,2),xyz(id1,1),'ok','markerfacecolor','k');
        hold(axes4,'on');
        plot(axes4, xyz(id1,1),-xyz(id1,3),'sk','markerfacecolor','w');
        plot(axes4, hat(:,2),hat(:,1),'r','linewidth',1.5);
        plot(axes4, hat(:,1),-hat(:,3),'r','linewidth',1.5);
        hold(axes4, 'off');
    case 'V vs E'
        plot(axes3, xyz(:,2),xyz(:,1),'ok-','markerfacecolor','k');
        hold(axes3,'on');
        plot(axes3, xyz(:,2),-xyz(:,3),'sk-','markerfacecolor','w')
        hold(axes3,'off');   

        plot(axes4, xyz(id1,2),xyz(id1,1),'ok','markerfacecolor','k');
        hold(axes4,'on');
        plot(axes4, xyz(id1,2),-xyz(id1,3),'sk','markerfacecolor','w');
        plot(axes4, hat(:,2),hat(:,1),'r','linewidth',1.5);
        plot(axes4, hat(:,2),-hat(:,3),'r','linewidth',1.5);
        hold(axes4,'off');
     case 'V vs H'
        plot(axes3, xyz(:,2),xyz(:,1),'ok-','markerfacecolor','k');
        hold(axes3,'on');
        plot(axes3, sqrt(xyz(:,1).^2 + xyz(:,2).^2),-xyz(:,3),'sk-','markerfacecolor','w')
        hold(axes3,'off');    

        plot(axes4, xyz(id1,2),xyz(id1,1),'ok','markerfacecolor','k');
        hold(axes4,'on');
        plot(axes4, sqrt(xyz(id1,1).^2 + xyz(id1,2).^2),-xyz(id1,3),'sk','markerfacecolor','w');
        plot(axes4, hat(:,2),hat(:,1),'r','linewidth',1.5);
        plot(axes4, sqrt(hat(:,1).^2 + hat(:,2).^2),-hat(:,3),'r','linewidth',1.5);
        hold(axes4,'off');
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
run_call = false;
list_call = true;
set(hObject,'UserData',list_call);
set(handles.Run,'UserData',run_call);

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
xin = [x, y, Dec, Inc];
id0=find(xin(:,1)<16);
id1=find(xin(:,1)>16);
xyz=ID2XYZ(deg2rad(xin(:,4)),deg2rad(xin(:,3)));
xyz=bsxfun(@times,xyz,xin(:,2));
tmp = xyz(id1,:)
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
[u,mad] = pmag_pca(x,y,z);
[coeff,score]=princomp(xyz(id1,:));

hat=[min(score(:,1));max(score(:,1))]*coeff(:,1)';
hat=bsxfun(@plus,hat,mean(xyz(id1,:)));

%%project data from different directions
plot_pca(xyz,hat, handles.axes3, handles.axes4, handles.popup_direction,id1);

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

xin=[0.00 9.2830e-05 339.9 57.9
2.50 7.5820e-05 325.7 49.1
5.00 6.2920e-05 321.3 45.9
10.00 5.2090e-05 314.8 41.7
15.00 4.4550e-05 310.3 38.7
20.00 3.9540e-05 305.0 37.0
30.00 3.2570e-05 303.9 34.7
40.00 2.5670e-05 303.0 32.3
50.00 2.2520e-05 303.6 32.4
60.00 1.9820e-05 299.8 30.8
70.00 1.3890e-05 292.5 31.0
80.00 1.2570e-05 297.0 25.6
90.00 0.5030e-05 299.3 11.3];
x = xin(:,1);
y = xin(:,2);
Dec = xin(:,3);
Inc = xin(:,4);

%x,y here show the values of first diagram.
% x = rawdata1.AFD0x28mT0x29;
% y = rawdata1.Intensity0x28tray0x29;

%x_pos, y_pos, z_pos here show the position in 3d space.
% x_pos = rawdata1.XtrayC;
% y_pos = rawdata1.YtrayC;
% z_pos = rawdata1.ZtrayC;

% Dec = rawdata1.Declination0x28tray0x29;
% Inc = rawdata1.Inclination0x28tray0x29;

% Dat = struct('name', name, 'x', x, 'y', y, 'x_pos', x_pos, 'y_pos', y_pos, 'z_pos', z_pos, 'Inc', Inc, 'Dec', Dec);
Dat = struct('name', name, 'x', x, 'y', y,'Inc', Inc, 'Dec', Dec);

set(hObject, 'UserData', Dat);
n = size(Dec);
set(handles.listbox1,'String', x);


% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


run_call = true;
list_call = false;
set(hObject,'UserData',run_call);
set(handles.listbox1,'UserData',list_call);

%%get the data info from load_dataset
dat_h = findobj('Tag', 'Load_dataset');
dat = get(dat_h,'UserData');

%%plot the first diagram
ax = findobj('Tag', 'axes1');
plot_first(dat.x,dat.y, handles.axes1);

%%plot the second diagram
plot_second(dat, handles.axes2);

%%plot the third diagram
xin = [dat.x, dat.y, dat.Dec, dat.Inc];
id0=find(xin(:,1)<16);
id1=find(xin(:,1)>16);
xyz=ID2XYZ(deg2rad(xin(:,4)),deg2rad(xin(:,3)));
xyz=bsxfun(@times,xyz,xin(:,2));
tmp = xyz(id1,:)
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);

[coeff,score]=princomp(xyz(id1,:));

hat=[min(score(:,1));max(score(:,1))]*coeff(:,1)';
hat=bsxfun(@plus,hat,mean(xyz(id1,:)));

%%project data from different directions
plot_pca(xyz,hat, handles.axes3, handles.axes4, handles.popup_direction,id1);



% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listbox1, 'String', '');
cla(handles.axes1,'reset');
cla(handles.axes2,'reset');
cla(handles.axes3,'reset');



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
run_h = findobj('Tag','Run');
run_call = get(run_h,'UserData');
list_h = findobj('Tag','listbox1');
list_call = get(list_h,'UserData');
if list_call
    data_list = getappdata(handles.listbox1,'data');
    plot_first(data_list.x, data_list.y, handles.axes1);
    
    plot_second(data_list, handles.axes2);
    
    xin = [data_list.x, data_list.y, data_list.Dec, data_list.Inc];
    id0=find(xin(:,1)<16);
    id1=find(xin(:,1)>16);
    xyz=ID2XYZ(deg2rad(xin(:,4)),deg2rad(xin(:,3)));
    xyz=bsxfun(@times,xyz,xin(:,2));
    tmp = xyz(id1,:)
    x = tmp(:,1);
    y = tmp(:,2);
    z = tmp(:,3);
    [coeff,score]=princomp(xyz(id1,:));
    hat=[min(score(:,1));max(score(:,1))]*coeff(:,1)';
    hat=bsxfun(@plus,hat,mean(xyz(id1,:)));
    plot_pca(xyz,hat, handles.axes3, handles.axes4, handles.popup_direction,id1);

end

if run_call  
    dat_h = findobj('Tag', 'Load_dataset');
    data_run = get(dat_h,'UserData');
    plot_first(data_run.x, data_run.y, handles.axes1);
    plot_second(data_run, handles.axes2);
    
    xin = [data_run.x, data_run.y, data_run.Dec, data_run.Inc];
    id0=find(xin(:,1)<16);
    id1=find(xin(:,1)>16);
    xyz=ID2XYZ(deg2rad(xin(:,4)),deg2rad(xin(:,3)));
    xyz=bsxfun(@times,xyz,xin(:,2));
    tmp = xyz(id1,:)
    x = tmp(:,1);
    y = tmp(:,2);
    z = tmp(:,3);
    [coeff,score]=princomp(xyz(id1,:));
    hat=[min(score(:,1));max(score(:,1))]*coeff(:,1)';
    hat=bsxfun(@plus,hat,mean(xyz(id1,:)));
    plot_pca(xyz,hat, handles.axes3, handles.axes4, handles.popup_direction,id1);

end


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
