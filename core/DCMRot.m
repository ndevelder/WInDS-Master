function [y,A]=DCMRot(x,t,A,rotseq,rev)
%% [y,A]=DCMRot(x,t,rotseq) -> Vector Rotation.
%
% Function performs a series of rotations about user-defined axes by
% user-defined angles over a series of vectors.
%
% ****Input(s)****
% x        1x3 (or Nx3) vector (array of vectors) to be rotated
% t        NxM array of rotation angles, where M=1..M corresponds to 
%          1st-Mth rotation order (degrees)
% A        Nx9 array representing preceeding rotation matrix
% rotseq   String (length M)indicating order of rotation sequence (Example:
%          'xyzy' indicates a rotation first about the x-axis, then y, then
%          z, then y
% rev      Compute transpose of DCM, then compute reverse sequence (if=1)
% 
% ****Output(s)****
% y        Nx3 array of rotated vectors
% A        Nx9 array representing rotation matrix
% 
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited February 26, 2010
%

%% Generate direction cosine matrix for rotation sequence
if isempty(A)
    A=zeros(size(t,1),9); %Form an identity array
    A(:,1:4:9)=1;
end

%Generate diagonal 1's and off-diagonal 0's
f0=zeros(size(t,1),1);
f1=ones(size(t,1),1);

%Speed up calculations by computing trig functions once
sint=sind(t);
cost=cosd(t);

for c1=1:length(rotseq) %Loop over number of rotation sequences
    if strcmpi(rotseq(c1),'x')
        R=[f1 f0 f0 f0 cost(:,c1) -sint(:,c1) f0 sint(:,c1) cost(:,c1)];
    elseif strcmpi(rotseq(c1),'y')
        R=[cost(:,c1) f0 sint(:,c1) f0 f1 f0 -sint(:,c1) f0 cost(:,c1)];
    elseif strcmpi(rotseq(c1),'z')
        R=[cost(:,c1) -sint(:,c1) f0 sint(:,c1) cost(:,c1) f0 f0 f0 f1]; 
    end
   
    B(:,1)=sum(R(:,1:3).*A(:,1:3:7),2);
    B(:,2)=sum(R(:,1:3).*A(:,2:3:8),2);
    B(:,3)=sum(R(:,1:3).*A(:,3:3:9),2);
    B(:,4)=sum(R(:,4:6).*A(:,1:3:7),2);
    B(:,5)=sum(R(:,4:6).*A(:,2:3:8),2);
    B(:,6)=sum(R(:,4:6).*A(:,3:3:9),2);
    B(:,7)=sum(R(:,7:9).*A(:,1:3:7),2);
    B(:,8)=sum(R(:,7:9).*A(:,2:3:8),2);
    B(:,9)=sum(R(:,7:9).*A(:,3:3:9),2);
    A=B;
end

if rev==1 %Compute transpose of DCM to reverse rotation sequence
    B(:,1)=A(:,1);
    B(:,2)=A(:,4);
    B(:,3)=A(:,7);
    B(:,4)=A(:,2);
    B(:,5)=A(:,5);
    B(:,6)=A(:,8);
    B(:,7)=A(:,3);
    B(:,8)=A(:,6);
    B(:,9)=A(:,9);
    A=B;
end

%% Apply rotation sequence to vector elements
if size(x,1)<size(A,1) %If a single vector undergoing a series of rotation, expand for 
                       %index multiplication
    x=repmat(x,size(A,1),1);
elseif size(x,1)>size(A,1) %If a single rotation seq. applied to multiple vectors, expand 
                           %for index multiplication
    A=repmat(A,size(x,1),1);
end

y(:,1)=sum(A(:,1:3).*x(:,1:3),2);
y(:,2)=sum(A(:,4:6).*x(:,1:3),2);
y(:,3)=sum(A(:,7:9).*x(:,1:3),2);