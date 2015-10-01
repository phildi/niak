function mat = niak_quat2mat(quat)
% Convert the quaternion info of a NIFTI header into an affine transformation
%
% SYNTAX:
% 
% MAT = NIAK_QUAT2MAT(QUAT)
% _________________________________________________________________________
% INPUT:
%
% QUAT
%   (structure) usually the HDR.DETAILS structure generated by NIAK_READ_VOL
%   with nifti (or analyze) file formats. The following fields are expected:
%
%   QUATERN_{B,C,D}
%      (scalar) the b, c and d parameters of the quaternion. 
%
%   QOFFSET_{X,Y,Z} 
%      (scalar) the translation 
%
%   PIXDIM
%      (vector) entries 2 to 4 define the size of the voxels.
%  
% _________________________________________________________________________
% OUTPUTS:
%
% ROT 
%       (array 3*N) the rotation parameters (in x, y and z planes). 
%       Unit is degrees.
%
% TSL 
%       (array 3*N) the translation parameters.
%
% _________________________________________________________________________
% COMMENTS:
%
% See the following web pages for more details
% http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/quatern.html
% http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html
%
% Pierre Bellec, 
% Research Centre of the Montreal Geriatric Institute
% & Department of Computer Science and Operations Research
% University of Montreal, Qubec, Canada, 2012
% Maintainer : pierre.bellec@criugm.qc.ca
% See licensing information in the code.
% Keywords : quaternion, nifti
%
% See licensing information in the code

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

b = quat.quatern_b;
c = quat.quatern_c;
d = quat.quatern_d;
a = sqrt(1-b^2-c^2-d^2);
ox = quat.qoffset_x;
oy = quat.qoffset_y;
oz = quat.qoffset_z;
sc = quat.pixdim(2:4);
if isfield(quat,'qfac')
    sc(end) = sc(end)*quat.qfac;
end
mat = zeros(3,4);
mat(4,4) = 1;
mat(1:3,1:3) = [ a*a+b*b-c*c-d*d , 2*b*c-2*a*d     ,  2*b*d+2*a*c     ; ...
                 2*b*c+2*a*d     , a*a+c*c-b*b-d*d ,  2*c*d-2*a*b     ; ...
                 2*b*d-2*a*c     , 2*c*d+2*a*b     ,  a*a+d*d-c*c-b*b ];
mat(1:3,1:3) = mat(1:3,1:3)*diag(sc);
mat(1:3,4) = [ox ; oy ; oz];
