function [img, head]=readseriesDCM(imnum,pathn)
% READSERIESDCM reads in a series of GE DICOM MRI images taken directly from the scanner.
% On a typical scanner, these images reside in the following directories:
%   /export/home1/sdc_image_pool/images
%
% [MR_IMAGE, MR_HEADER]=READSERIESDCM([first image last image], path)
% [MR_IMAGE, MR_HEADER]=READSERIESDCM([image_list], path)
% 
% Function to read in images from an MR image series in 
% GEMS DICOM format.
% Images saved in a 3D variable (img).  Header information 
% stored in 3D variable (head).
%
% Optional arguments are first image and last image which
% specify the start/stop parameters for reading images.
% Each argument must be an integer.  If these parameters are not
% specified, the ENTIRE series will be read.
%
% The path is also an optional argument.  If specified, the user
% interface selection of a study is bypassed.
%
% SEE ALSO:  READLX, READHEADERLX, READSIGNA, READSERIES, READHEADER, READOFFSET, DBPATH, INAME

% Author: R. Jason Stafford
% Date: 01/05
% Revision:	(based on READSERIES and READSERIESLX)
%		

% If the images aren't specified, allow user to locate

if nargin<2,[pathn filen] = dbpath;,end

if size(dir([pathn '*.dcm']),1)>0 % Nordic-exported images: 1.dcm, 2.dcm, etc.
   tmp = dir([pathn '*.dcm']);    
    fnames = strvcat(tmp.name);
    nimgs = size(fnames,1);

    % Resorts fnames
    for qq = 1:nimgs
        c = strfind(fnames(qq,:),'.')-1; % c = index where file #'s end (where the dot is)
        sortnumbers(qq,1) = str2num(fnames(qq,(1:c))); % extracts starting #s
    end
    [sortnumbers I] = sort(sortnumbers); % how they should be re-sorted
    fnames = fnames(I,:); % Resort them

elseif size(dir([pathn '*.dcm.*']),1)>0 % Nordic-exported images: 1.dcm, 2.dcm, etc.
   tmp = dir([pathn '*.dcm.*']);    
    fnames = strvcat(tmp.name);
    nimgs = size(fnames,1);
    % headtmp = readheaderlx(pathn,fnames(1,:));

    % Resorts fnames
    for qq = 1:nimgs
        c = strfind(fnames(qq,:),'.')-1; % c = index where file #'s end (where the dot is)
        sortnumbers(qq,1) = str2num(fnames(qq,(1:c))); % extracts starting #s
    end
    [sortnumbers I] = sort(sortnumbers); % how they should be re-sorted
    fnames = fnames(I,:); % Resort them

elseif size(dir([pathn '*IMG*']),1)>0 % Olea-exported images: IMG0, IMG1, etc.
    tmp=dir([pathn '*IMG*']);, 
    fnames = strvcat(tmp.name);
    nimgs = size(fnames,1);
    % headtmp = readheaderlx(pathn,fnames(1,:));

    % Resorts fnames
    for qq = 1:nimgs
        sortnumbers(qq,1) = str2num(fnames(qq,(4:end))); % extracts ending #s
    end
    [sortnumbers I] = sort(sortnumbers); % how they should be re-sorted
    fnames = fnames(I,:); % Resort them
    
elseif size(dir([pathn '*IM*']),1)>0 % AW-exported images: IM0, IM1, etc.
    tmp=dir([pathn '*IM*']);, 
    fnames = strvcat(tmp.name);
    nimgs = size(fnames,1);
    % headtmp = readheaderlx(pathn,fnames(1,:));

    % Resorts fnames
    for qq = 1:nimgs
        sortnumbers(qq,1) = str2num(fnames(qq,(3:end))); % extracts ending #s
    end
    [sortnumbers I] = sort(sortnumbers); % how they should be re-sorted
    fnames = fnames(I,:); % Resort them
    
elseif size(dir([pathn '*.MRDC.*']),1)>0  % GE Scanner images
   tmp = dir([pathn '*.MRDC.*']);
    fnames = strvcat(tmp.name);
    nimgs = size(fnames,1);
    % headtmp = readheaderlx(pathn,fnames(1,:));

    % Resorts fnames
    for qq = 1:nimgs
        c = strfind(fnames(qq,:),'.MRDC.')+6; % c = index where file #'s start 
        sortnumbers(qq,1) = str2num(fnames(qq,(c:end))); % extracts #s
    end
    [sortnumbers I] = sort(sortnumbers); % how they should be re-sorted
    fnames = fnames(I,:); % Resort them
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Read in images from a series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CASE 1:  If the image range isn't specified

if nargin == 0,         
    img1=1;
   imgf=size(fnames,1);
   if nargout>1,
       [img, head]=readseriesDCM([img1 imgf],pathn);
   else
       img=readseriesDCM([img1 imgf],pathn);
   end

% CASE 2:  Function is called with initial and final image numbers
else
    %disp([' '])
   % disp(['Images found:'])
   % disp([' '])
    % disp([strvcat(tmp.name)])
    %disp([' '])
    %disp(['Reading images ...'])
    disp([' '])
   if length(imnum)==2,
      imnum = [imnum(1):imnum(2)];
   end
   
   jj=1;
   for i = imnum,
      [tmpi, tmph]=readDCM(pathn,deblank(fnames(i,:)));
      img(:,:,jj)=tmpi;
      if nargout>1,head(jj)=tmph;,end
      jj = jj+1;
	end
end

% end function