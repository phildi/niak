% This is a script to demonstrate the usage of :
% NIAK_BRICK_MOTION_CORRECTION
%
% SYNTAX:
% Just type in NIAK_DEMO_MOTION_CORRECTION 
%
% OUTPUT:
%
% This script will clear the workspace !!
% It will apply a motion correction on the functional data of subject
% 1 and use the default output name.
%
% Note that the path to access the demo data is stored in a variable
% called GB_NIAK_PATH_DEMO defined in the NIAK_GB_VARS script.
% 
% The demo database exists in multiple file formats. By default, it is
% using 'minc2' files. You can change that by changing the variable
% GB_NIAK_FORMAT_DEMO in the file NIAK_GB_VARS.
%
% Copyright (c) Pierre Bellec, Montreal Neurological Institute, 2008.
% Maintainer : pbellec@bic.mni.mcgill.ca
% See licensing information in the code.
% Keywords : medical imaging, slice timing, fMRI

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

clear
niak_gb_vars

%% Setting input/output files
switch gb_niak_format_demo
        
    case 'minc2' % If data are in minc2 format
        
        %% The two datasets have actually been acquired in the same
        %% session, but this is just to demonstrate how the procedure works
        %% in general.
        files_in.fmri = cat(2,gb_niak_path_demo,filesep,'func_motor_subject2.mnc');
        files_in.target = cat(2,gb_niak_path_demo,filesep,'func_motor_subject2_mean.mnc');
    case 'minc1'
        
        %% The two datasets have actually been acquired in the same
        %% session, but this is just to demonstrate how the procedure works
        %% in general.
        files_in.fmri = cat(2,gb_niak_path_demo,filesep,'func_motor_subject2.mnc.gz');
        files_in.target = cat(2,gb_niak_path_demo,filesep,'func_motor_subject2_mean.mnc.gz');
        
    otherwise 
        
        error('niak:demo','%s is an unsupported file format for this demo. See help to change that.',gb_niak_format_demo)        
end

%% Setting output files
files_out = ''; % use default names

%% Options

opt.flag_test = 0; % Actually perform the motion correction
opt.tol = 0.0005;
opt.fwhm = 5;

%% Generate the target volume
[hdr,vol] = niak_read_vol(files_in.fmri);
hdr.file_name = files_in.target;
niak_write_vol(hdr,vol(:,:,:,1));

[files_in,files_out,opt] = niak_brick_motion_correction_dev(files_in,files_out,opt);

tab = niak_read_tab(files_out);
subplot(1,2,2)
plot(tab(:,1:3));
subplot(1,2,1)
plot(tab(:,4:6));

%% Note that opt.interpolation_method has been updated, as well as files_out

