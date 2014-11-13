function [mrimg, header] = readDCM(path,file)
% READDCM reads in a single GE EXCITE MRI image taken directly 
% from the scanner (/export/home1/sdc_image_pool/images) running
% on the Excite platform (GEDM images).
% On a typical scanner, these images reside in the following directories:
%				(fill in later)
% [MR_IMAGE, MR_HEADER]=READDCM(path,filename)
% [MR_IMAGE, MR_HEADER]=READDCM;
% 
% Function to read in images from an MR image series in 
% GEMS Lightning format.
% Images saved in a 3D variable (img).  Header information 
% stored in 3D variable (head).
%
% The path and filename are also an optional arguments.  If specified, the user
% interface selection of a study is bypassed.
%
% SEE ALSO:  READLX, READHEADERLX, READSIGNA, READSERIES, READHEADER, READOFFSET, DBPATH, INAME

% Author: R. Jason Stafford
% Date: 01/01
% Revision:	(based on READSIGNA)

% Share some info
% global header
% global mrimg

if nargin==0, [path, file]=dbpath; end

% Read in header info
header=dicominfo([path file],'dictionary','gems-dicom-dict.txt');
%header=dicominfo([path file]);

   
% Read in the GE DICOM
mrimg=dicomread([path file]);

disp(['Read ' path file])