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

% Last Modified by GUIDE v2.5 26-Feb-2016 09:33:24

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
handles.laxis = axes('parent',hObject,'units','normalized','position',[0 0 1 1],'visible','off');
lbls = findobj(hObject,'-regexp','tag','latex_*');
display(length(lbls));
for i=1:length(lbls)
      l = lbls(i);
      % Get current text, position and tag
      set(l,'units','normalized');
      s = get(l,'string');
      p = get(l,'position');
      t = get(l,'tag');
      % Remove the UICONTROL
      delete(l);
      % Replace it with a TEXT object 
      handles.(t) = text(p(1),p(2),s,'interpreter','latex');
end

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


function plot_bpca(xyz,MAP_I,MAP_D,VN,VE,VH,Emu,popup_direction,axes4)
MAP_I= deg2rad(MAP_I);
MAP_D = deg2rad(MAP_D);
map_xyz = ID2XYZ(MAP_I, MAP_D);

min_pt = min(VN.H{1}(:,2))*map_xyz;
max_pt = max(VN.H{1}(:,2))*map_xyz;
hat=bsxfun(@plus,[min_pt;max_pt],Emu');

%scale the VN.H line
siz = size(VN.H{1});
len = siz(1);

k1 = (hat(1,1)-hat(2,1))/(hat(1,2)-hat(2,2));
b1 = hat(1,1)-k1*hat(1,2);
x1_min = ((VN.H{1}(1,2)+VN.H{2}(1,2))/2 - b1)/k1;
x1_max = ((VN.H{1}(len,2)+VN.H{2}(len,2))/2 - b1)/k1;

%scale the VN.V line
k2 = (-hat(1,3)+hat(2,3))/(hat(1,1)-hat(2,1));
b2 = -hat(1,3)-k2*hat(1,1);
x2_min = ((VN.V{1}(1,2)+VN.V{2}(1,2))/2 - b2)/k2;
x2_max = ((VN.V{1}(len,2)+VN.V{2}(len,2))/2 - b2)/k2;

%scale VE.H line
k3 = k1;
b3 = b1;
x3_min = ((VE.H{1}(1,2)+VE.H{2}(1,2))/2 - b3)/k3;
x3_max = ((VE.H{1}(len,2)+VE.H{2}(len,2))/2 - b3)/k3;

%scale VE.V line
k4 = (-hat(1,3)+hat(2,3))/(hat(1,2)-hat(2,2));
b4 = -hat(1,3) - k4*hat(1,2);
x4_min = ((VE.V{1}(1,2)+VE.V{2}(1,2))/2 - b4)/k4;
x4_max = ((VE.V{1}(len,2)+VE.V{2}(len,2))/2 - b4)/k4;

%scale VH.H line
k5 = k1;
b5 = b1;
x5_min = ((VH.H{1}(1,2)+VH.H{2}(1,2))/2 - b5)/k5;
x5_max = ((VH.H{1}(len,2)+VH.H{2}(len,2))/2 - b5)/k5;

%scale VH.V line
k6 = (-hat(1,3)+hat(2,3))/(sqrt(hat(1,1).^2+hat(1,2).^2) - sqrt(hat(2,1).^2+hat(2,2).^2));
b6 = -hat(1,3)-k6*sqrt(hat(1,1).^2+hat(1,2).^2);
x6_min = ((VH.V{1}(1,2)+VH.V{2}(1,2))/2 - b6)/k6;
x6_max = ((VH.V{1}(len,2)+VH.V{2}(len,2))/2 - b6)/k6;



contents = get(popup_direction,'String');
content = contents{get(popup_direction,'Value')};
switch content
    case 'V vs N'
        plot(axes4, VN.H{1}(:,1), VN.H{1}(:,2),'--k');
        hold(axes4,'on');
        plot(axes4,[x1_min;x1_max],[(VN.H{1}(1,2)+VN.H{2}(1,2))/2;(VN.H{1}(len,2)+VN.H{2}(len,2))/2]);
        plot(axes4,VN.H{2}(:,1), VN.H{2}(:,2),'--k');
        p1 = plot(axes4,xyz(:,2),xyz(:,1),'ok','markerfacecolor','k');

        plot(axes4,VN.V{1}(:,1), VN.V{1}(:,2),'--r');
        plot(axes4,VN.V{2}(:,1), VN.V{2}(:,2),'--r');
        plot(axes4,[x2_min;x2_max],[(VN.V{1}(1,2)+VN.V{2}(1,2))/2;(VN.V{1}(len,2)+VN.V{2}(len,2))/2]);
        p2 = plot(axes4,xyz(:,1),-xyz(:,3),'sk','markerfacecolor','w');
        legend([p2,p1],'Vertical','Horizontal','location','northeastoutside');
        set(axes4,'Dataaspectratio', [1 1 1]);
        plot_coord(axes4);
        hold(axes4,'off');
    case 'V vs E'
        plot(axes4, VE.H{1}(:,1), VE.H{1}(:,2),'--k');
        hold(axes4,'on');
        plot(axes4,VE.H{2}(:,1),VE.H{2}(:,2),'--k');
        plot(axes4,[x3_min;x3_max],[(VE.H{1}(1,2)+VE.H{2}(1,2))/2;(VE.H{1}(len,2)+VE.H{2}(len,2))/2]);
        p1 = plot(axes4,xyz(:,2),xyz(:,1),'ok','markerfacecolor','k');

        plot(axes4,VE.V{1}(:,1), VE.V{1}(:,2),'--r');
        plot(axes4,VE.V{2}(:,1), VE.V{2}(:,2),'--r');
        plot(axes4,[x4_min;x4_max],[(VE.V{1}(1,2)+VE.V{2}(1,2))/2;(VE.V{1}(len,2)+VE.V{2}(len,2))/2]);
        p2 = plot(axes4,xyz(:,2),-xyz(:,3),'sk','markerfacecolor','w');
        legend([p2,p1],'Vertical','Horizontal','location','northeastoutside');
        set(axes4,'Dataaspectratio', [1 1 1]);
        plot_coord(axes4);
        hold(axes4,'off');        
    case 'V vs H'
        plot(axes4, VH.H{1}(:,1), VH.H{1}(:,2),'--k');
        hold(axes4,'on');
        plot(axes4,VH.H{2}(:,1), VH.H{2}(:,2),'--k');
        plot(axes4,[x5_min;x5_max],[(VH.H{1}(1,2)+VH.H{2}(1,2))/2;(VH.H{1}(len,2)+VH.H{2}(len,2))/2]);
        p1 = plot(axes4,xyz(:,2),xyz(:,1),'ok','markerfacecolor','k');
        % 
        plot(axes4,VH.V{1}(:,1), VH.V{1}(:,2),'--r');
        plot(axes4,VH.V{2}(:,1), VH.V{2}(:,2),'--r');
        plot(axes4,[x6_min;x6_max],[(VH.V{1}(1,2)+VH.V{2}(1,2))/2;(VH.V{1}(len,2)+VH.V{2}(len,2))/2]);
        p2 = plot(axes4, sqrt(xyz(:,1).^2 + xyz(:,2).^2),-xyz(:,3),'sk','markerfacecolor','w');
        legend([p2,p1],'Vertical','Horizontal','location','northeastoutside');
        set(axes4,'Dataaspectratio', [1 1 1]);
        plot_coord(axes4);
        hold(axes4,'off');

end

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
        p1 = plot(axes4, xyz(:,2),xyz(:,1),'ok','markerfacecolor','k');
        hold(axes4,'on');
        p2 = plot(axes4, xyz(:,1),-xyz(:,3),'sk','markerfacecolor','w');
        legend([p2,p1],'Vertical','Horizontal','location','northeastoutside');
        plot(axes4, hat(:,2),hat(:,1),'r','linewidth',1.5);
        plot(axes4, hat(:,1),-hat(:,3),'r','linewidth',1.5);
        set(axes4,'Dataaspectratio', [1 1 1]);
        plot_coord(axes4);
        hold(axes4, 'off');
    case 'V vs E'
        p1 = plot(axes4, xyz(:,2),xyz(:,1),'ok','markerfacecolor','k');
        hold(axes4,'on');
        p2 = plot(axes4, xyz(:,2),-xyz(:,3),'sk','markerfacecolor','w');
        legend([p2,p1],'Vertical','Horizontal','location','northeastoutside');
        plot(axes4, hat(:,2),hat(:,1),'r','linewidth',1.5);
        plot(axes4, hat(:,2),-hat(:,3),'r','linewidth',1.5);
        set(axes4,'Dataaspectratio', [1 1 1]);
        plot_coord(axes4);
        hold(axes4,'off');
     case 'V vs H'
        p1 = plot(axes4, xyz(:,2),xyz(:,1),'ok','markerfacecolor','k');
        hold(axes4,'on');
        p2 = plot(axes4, sqrt(xyz(:,1).^2 + xyz(:,2).^2),-xyz(:,3),'sk','markerfacecolor','w');
        legend([p2,p1],'Vertical','Horizontal','location','northeastoutside');
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
        p1 = plot(axes3, xyz(:,2),xyz(:,1),'ok-','markerfacecolor','k');
        hold(axes3,'on');
        p2 = plot(axes3, xyz(:,1),-xyz(:,3),'sk-','markerfacecolor','w');
        legend([p2,p1],'Vertical','Horizontal','location','northeastoutside');
        plot(axes3, color_xyz(:,2),color_xyz(:,1),'*', 'color', 'r');
        plot(axes3, color_xyz(:,1),-color_xyz(:,3),'*', 'color', 'b');
        
%         set(axes3,'Dataaspectratio', [1 1 1]);
        plot_coord(axes3);
        hold(axes3,'off');
    case 'V vs E'
        p1 = plot(axes3, xyz(:,2),xyz(:,1),'ok-','markerfacecolor','k');
        hold(axes3,'on');
        p2 = plot(axes3, xyz(:,2),-xyz(:,3),'sk-','markerfacecolor','w');
        legend([p2,p1],'Vertical','Horizontal','location','northeastoutside');
        plot(axes3, color_xyz(:,2),color_xyz(:,1),'*', 'color', 'r');
        plot(axes3, color_xyz(:,2),-color_xyz(:,3),'*', 'color', 'b');
        plot_coord(axes3);
        set(axes3,'Dataaspectratio', [1 1 1]);
        hold(axes3,'off');   
     case 'V vs H'
        p1 = plot(axes3, xyz(:,2),xyz(:,1),'ok-','markerfacecolor','k');
        hold(axes3,'on');
        p2 = plot(axes3, sqrt(xyz(:,1).^2 + xyz(:,2).^2),-xyz(:,3),'sk-','markerfacecolor','w');
        legend([p2,p1],'Vertical','Horizontal','location','northeastoutside');
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
    if i == 1
        plot(axes2, [-0.5;0.5],[-sqrt(3)/2;sqrt(3)/2],'k');
        plot(axes2, [-sqrt(3)/2;sqrt(3)/2],[-0.5;0.5],'k');
        plot(axes2, [-1;1],[0;0],'k');
        plot(axes2, [-sqrt(3)/2;sqrt(3)/2],[0.5;-0.5],'k');
        plot(axes2, [-0.5;0.5],[sqrt(3)/2;-sqrt(3)/2],'k');
        plot(axes2, [0;0],[-1;1],'k');
    end
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
setappdata(handles.Load_dataset,'bayesian',false); %Load_dataset.bayesian to decide whether BPCA is being used.
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

set(handles.MAD, 'String', floor(mad*100)/100);
[I,D] = XYZ2ID(u');
set(handles.Inclination, 'String' ,floor(rad2deg(I)*100)/100); %only display 2 digits on the right of decimal
set(handles.Declination, 'String', ceil(rad2deg(D)*100)/100);


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
  'Select Data File');
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

x = rawdata1.AFD0x28mT0x29;
y = rawdata1.Intensity0x28tray0x29;


Dec = rawdata1.Declination0x28tray0x29;
Inc = rawdata1.Inclination0x28tray0x29;

Dat = struct('name', name, 'x', x, 'y', y,'Inc', Inc, 'Dec', Dec);

set(hObject, 'UserData', Dat); %this Load_dataset.userdata will be the data input
n = size(Dec);
set(handles.listbox1,'String', x);
set(handles.sampleID,'String',name(1,:));
setappdata(hObject,'sampleID',name(1,:));

% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%reset everything
set(handles.listbox1, 'String', '');
cla(handles.axes1,'reset');
cla(handles.axes2,'reset');
cla(handles.axes3,'reset');
cla(handles.axes4,'reset');
set(handles.Load_dataset, 'UserData', []);
setappdata(handles.listbox1, 'data', []);
setappdata(hObject, 'bayesian',false);

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

data_list = getappdata(handles.listbox1,'data'); % get the selected data from listbox1
plot_first(data_list.x, data_list.y, handles.axes1);

plot_second(data_list, handles.axes2);
dat_h = findobj('Tag', 'Load_dataset');
data_run = get(dat_h,'UserData');     % get the data input from Load_dataset
plot_third(data_run, data_list, handles.axes3, handles.popup_direction);
bayesian = getappdata(handles.Load_dataset, 'bayesian');
if bayesian == false    
    plot_pca(data_list, handles.axes4, handles.popup_direction);
else
    % get the prameter for plotting BPCA diagram from bayesian, rather than go through BPCA algorithm again
    MAP_I = getappdata(handles.Bayesian,'MAP_I');
    MAP_D = getappdata(handles.Bayesian,'MAP_D');
    theta95 = getappdata(handles.Bayesian,'theta95');
    MD0 = getappdata(handles.Bayesian,'MD0');
    VN = getappdata(handles.Bayesian,'VN');
    VE = getappdata(handles.Bayesian,'VE');
    VH = getappdata(handles.Bayesian,'VH');
    Emu = getappdata(handles.Bayesian,'Emu');
    xyz = getappdata(handles.Bayesian,'xyz');
    plot_bpca(xyz,MAP_I,MAP_D,VN,VE,VH,Emu,handles.popup_direction,handles.axes4);
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
try
    x = data.x;
catch
    ed = errordlg('there must be at least two data point as input', 'Error');
    set(ed, 'WindowStyle', 'modal');
    uiwait(ed);    
end

setappdata(handles.Load_dataset,'bayesian',true);
xyz=ID2XYZ(deg2rad(data.Inc),deg2rad(data.Dec));
xyz=bsxfun(@times,xyz,data.y);
p = str2double(get(handles.confidence,'String'));
if p>=1 | p<=0 | isnan(p)
    ed = errordlg('confidence range is between 0 and 1', 'Error');
    set(ed, 'WindowStyle', 'modal');
    uiwait(ed);
else


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

plot_bpca(xyz,MAP_I,MAP_D,VN,VE,VH,Emu,handles.popup_direction,handles.axes4);
setappdata(hObject,'MAP_I',MAP_I);
setappdata(hObject,'MAP_D',MAP_D);
setappdata(hObject,'theta95',theta95);
setappdata(hObject,'MD0',MD0);
setappdata(hObject,'VN',VN);
setappdata(hObject,'VE',VE);
setappdata(hObject,'VH',VH);
setappdata(hObject,'Emu',Emu);
setappdata(hObject,'xyz',xyz);

%display parameters

set(handles.MAP_I,'String',floor(MAP_I*100)/100);
set(handles.MAP_D,'String',floor(MAP_D*100)/100);
set(handles.Ip1,'String',floor(Ip(1)*100)/100);
set(handles.Dp1,'String',floor(Dp(1)*100)/100);
set(handles.Ip2,'String',floor(Ip(2)*100)/100);
set(handles.Dp2,'String',floor(Dp(2)*100)/100);
set(handles.mD0,'String',floor(MD0*100)/100);
set(handles.theta,'String',floor(theta95*100)/100);

close(h);
set(handles.listbox1,'Enable','on');
set(handles.Load_dataset,'Enable','on');
set(handles.clear,'Enable','on');
set(handles.Bayesian,'Enable','on');
set(handles.popup_direction,'Enable','on');

end



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


% --- Executes on button press in output.
function output_Callback(hObject, eventdata, handles)
% hObject    handle to output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile('*.csv', 'Save output as');
setappdata(handles.append_to_output,'filename',filename);
setappdata(handles.append_to_output,'pathname',pathname);
header0 = 'SampleID';
header1 = 'PCA_I';
header2 = 'PCA_D';
header3 = 'MAD';
header4 = 'MAP_I';
header5 = 'MAP_D';
header6 = 'MD0';
header7 = 'theta';
header8 = 'I_p1';
header9 = 'I_p2';
header10 = 'D_p1';
header11 = 'D_p2';
header12 = 'confidence';
fid = fopen(strcat(pathname,filename),'w');
PCA_I = str2double(get(handles.Inclination,'String'));
PCA_D = str2double(get(handles.Declination,'String'));
MAD = str2double(get(handles.MAD,'String'));
sampleID = getappdata(handles.Load_dataset,'sampleID');
bayesian = getappdata(handles.Load_dataset,'bayesian');
if bayesian == false
    fprintf(fid, [header0 ',' header1 ',' header2 ',' header3 ',' header4 ',' header5 ',' header6 ',' header7 ',' header8 ',' header9 ',' header10 ',' header11 ',' header12 '\n']);
    fprintf(fid, '%s,%.2f,%.2f,%.2f,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', sampleID,PCA_I, PCA_D, MAD, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA','NA');
    fclose(fid);

else
    MAP_I = str2double(get(handles.MAP_I,'String'));
    MAP_D = str2double(get(handles.MAP_D,'String'));
    mD0 = str2double(get(handles.mD0,'String'));
    theta = str2double(get (handles.theta,'String'));
    Ip1 = str2double(get(handles.Ip1,'String'));
    Ip2 = str2double(get(handles.Ip2,'String'));
    Dp1 = str2double(get(handles.Dp1,'String'));
    Dp2 = str2double(get(handles.Dp2,'String'));
    confidence = str2double(get(handles.confidence,'String'));
    fprintf(fid, [header0 ',' header1 ',' header2 ',' header3 ',' header4 ',' header5 ',' header6 ',' header7 ',' header8 ',' header9 ',' header10 ',' header11 ',' header12 '\n']);
    fprintf(fid, '%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n', sampleID,PCA_I, PCA_D, MAD, MAP_I,MAP_D, mD0, theta, Ip1, Ip2, Dp1, Dp2,confidence);
    fclose(fid);
end

% --- Executes on button press in append_to_output.
function append_to_output_Callback(hObject, eventdata, handles)
% hObject    handle to append_to_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = getappdata(hObject,'filename');
pathname = getappdata(hObject,'pathname');
if ~exist(strcat(pathname,filename),'file')
    ed = errordlg('press Output to select output file first', 'Error');
    set(ed, 'WindowStyle', 'modal');
    uiwait(ed); 
end
    

fid = fopen(strcat(pathname,filename),'a');
PCA_I = str2double(get(handles.Inclination,'String'));
PCA_D = str2double(get(handles.Declination,'String'));
MAD = str2double(get(handles.MAD,'String'));
sampleID = getappdata(handles.Load_dataset,'sampleID');
bayesian = getappdata(handles.Load_dataset,'bayesian');
if bayesian == false
    fprintf(fid, '%s,%.2f,%.2f,%.2f,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', sampleID,PCA_I, PCA_D, MAD, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA','NA');
    fclose(fid);
else
    MAP_I = str2double(get(handles.MAP_I,'String'));
    MAP_D = str2double(get(handles.MAP_D,'String'));
    mD0 = str2double(get(handles.mD0,'String'));
    theta = str2double(get(handles.theta,'String'));
    Ip1 = str2double(get(handles.Ip1,'String'));
    Ip2 = str2double(get(handles.Ip2,'String'));
    Dp1 = str2double(get(handles.Dp1,'String'));
    Dp2 = str2double(get(handles.Dp2,'String'));
    confidence = str2double(get(handles.confidence,'String'));
    fprintf(fid, '%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n', sampleID,PCA_I, PCA_D, MAD, MAP_I,MAP_D, mD0, theta, Ip1, Ip2, Dp1, Dp2,confidence);
    fclose(fid);
end
