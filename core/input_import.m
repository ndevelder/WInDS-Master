function [airfoils,blade,turbine,platform,wind]=input_import(filename)
%% [airfoils,blade,turbine,platform,wind]=input_import(filename) -> FAST input files importer.
%
% Function imports FAST simulation input files
% 
% ****Input(s)****
% filename  String containing path to FAST input file (.fst)
% 
% ****Output(s)****
% airfoils  Structure containing airfoil performance tables
% blade     Structure containing blade geometry from FAST input file
% turbine   Structure containing turbine geometry from FAST input file
% platform  Structure containing platform geometry from FAST input file
% wind      Structure containing wind data file location
% 
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited February 23, 2010
%

%% Use FAST input file to ID other relevant files
fn=strread(char(filename),'%s','delimiter','\\');
fstfile=char(fn(end));
fstpath=filename(1:end-length(fstfile));
turbine.filename=[fstpath fstfile];
data=importdata(turbine.filename,'\t'); %Import FAST input file

% Identify platform property file
pf=sscanf(char(data(131)),'%i');
if pf>=2
platform.filename=strread(char(data(132)),'%s','delimiter','"');
platform.filename=[fstpath char(platform.filename(2))];
else
platform.filename='No platform model used.';
end

% Identify AeroDyn file
blade.filename=strread(char(data(161)),'%s','delimiter','"');
blade.filename=[fstpath char(blade.filename(2))];

%% Import turbine and blade properties from FAST input file
blade.TipRad=sscanf(char(data(78)),'%f');
blade.HubRad=sscanf(char(data(79)),'%f');
turbine.NumBl=sscanf(char(data(9)),'%i');
turbine.OverHang=sscanf(char(data(83)),'%f');
turbine.TowerHt=sscanf(char(data(87)),'%f');
turbine.Twr2Shft=sscanf(char(data(88)),'%f');
turbine.ShftTilt=sscanf(char(data(90)),'%f');
turbine.PreCone(1)=sscanf(char(data(92)),'%f');
turbine.PreCone(2)=sscanf(char(data(93)),'%f');
turbine.PreCone(3)=sscanf(char(data(94)),'%f');
clear data

%% Import AeroDyn file and individual airfoil files
% Identify TurbSim-based wind input file
data=importdata(blade.filename,'\t');
wind.filename=strread(char(data(10)),'%s','delimiter','"');
wind.filename=[fstpath char(wind.filename(2))];

% Count up number of airfoils and blade sections, import airfoil tables, 
% then import blade properties as a structure
nblades=sscanf(char(data(18)),'%i');
airfoils.Names=cell(nblades,1);
for c1=1:nblades
    af=strread(char(data(18+c1)),'%s','delimiter','"');
    af=char(af(2));
    adata=importdata([fstpath af],'\t');
    
    if(isfield(adata,'textdata')==0) %Sometimes will import a cell structure, check for this
        adata1=importdata([fstpath af],' ',14);
        adata2=importdata([fstpath af],'\t');
        adata2(15:end)=[];
        clear adata
        adata.data=adata1.data;
        adata.textdata=adata2;
    end

    id=isnan(adata.data(2,:));
    adata.data(:,id)=[];
    
    af=strread(char(af),'%s','delimiter','\\');
    af=char(af(end));
    
    airfoils.Names(c1,1)={genvarname(af(1:end-4))};
    
    id=find(diff(adata.data(:,1))==0); %ID non-distinct values for AoA
    adata.data(id,:)=[]; %#ok<FNDSB>
    
    eval(['airfoils.profiles(' num2str(c1) ' ,1).StallAoA=' ...
        'sscanf(char(adata.textdata(5)),''%f'');'])
    eval(['airfoils.profiles(' num2str(c1) ' ,1).Cn0AoA=' ...
        'sscanf(char(adata.textdata(9)),''%f'');'])
    eval(['airfoils.profiles(' num2str(c1) ' ,1).Lift0Cn=' ...
        'sscanf(char(adata.textdata(10)),''%f'');'])
    eval(['airfoils.profiles(' num2str(c1) ' ,1).StallAoACn=' ...
        'sscanf(char(adata.textdata(11)),''%f'');'])
    eval(['airfoils.profiles(' num2str(c1) ' ,1).StallAoANCn=' ...
        'sscanf(char(adata.textdata(12)),''%f'');'])
    eval(['airfoils.profiles(' num2str(c1) ' ,1).CdminAoA=' ...
        'sscanf(char(adata.textdata(13)),''%f'');'])
    eval(['airfoils.profiles(' num2str(c1) ' ,1).Cdmin=' ...
        'sscanf(char(adata.textdata(14)),''%f'');'])
    
    eval(['airfoils.profiles(' num2str(c1) ' ,1).AoA=adata.data(:,1);'])
    eval(['airfoils.profiles(' num2str(c1) ' ,1).Cl=adata.data(:,2);'])
    eval(['airfoils.profiles(' num2str(c1) ' ,1).Cd=adata.data(:,3);'])
    eval(['airfoils.profiles(' num2str(c1) ' ,1).Cm=adata.data(:,4);'])
    clear af adata adata1 adata2 id
end

dm=19+nblades;
ivnames=textscan(char(data(dm+1,:)),'%s');
ivnames=genvarname(cell(ivnames{1,1}));
data=char(data(dm+1:length(data),:));
ndata=zeros(size(data,1)-1,5);
for c2=2:size(data,1)
    ndata(c2-1,:)=sscanf(data(c2, :)', '%f %f %f %f %d', [1, inf]);
end
for c3=1:5
    eval(['blade.' char(ivnames(c3)) '=ndata(:,c3);'])
end

%% Import platform properties
if pf==0 || pf==1
    platform.Type='onshore';
    platform.TwrDraft=0;
    platform.PtfmCM=0;
    platform.PtfmRef=0;
    platform.PtfmDraft=0;
    platform.PtfmDiam=0;
elseif pf==2
    data=importdata(platform.filename,'\t');
    platform.Type='fixedoffshore';
    platform.TwrDraft=sscanf(char(data(19)),'%f');
    platform.PtfmCM=sscanf(char(data(20)),'%f');
    platform.PtfmRef=sscanf(char(data(21)),'%f');
    platform.PtfmDraft=sscanf(char(data(36)),'%f'); %Water depth
    platform.PtfmDiam=sscanf(char(data(31)),'%f');
elseif pf==3
    data=importdata(platform.filename,'\t');
    platform.Type=strread(char(data(29)),'%s','delimiter','"');
    platform.Type=strread(char(platform.Type(2)),'%s','delimiter','\\');
    platform.Type=char(platform.Type(end));
    platform.TwrDraft=sscanf(char(data(19)),'%f');
    platform.PtfmCM=sscanf(char(data(20)),'%f');
    platform.PtfmRef=sscanf(char(data(21)),'%f');
    platform.PtfmDraft=sscanf(char(data(32)),'%f');
    platform.PtfmDiam=sscanf(char(data(33)),'%f');
end