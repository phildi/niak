function [in,out,opt] = niak_brick_table_clusters(in,out,opt)
% Generate a .csv table summary of the spatial distribution of clusters.
%
% SYNTAX:
% [IN,OUT,OPT] = NIAK_BRICK_TABLE_CLUSTERS(IN,OUT,OPT)
%
% INPUTS:
%   IN.PART (string) a file name of a 3D partition (cluster I filled with Is).
%       Note that the cluster numerical values do not need to be 1, ... , N.
%       The background is filled with 0s.
%   IN.LABELS_PART (string, default 'clusterI') a .csv file name where row I 
%      has a (string) label for network I. The .csv can also have two columns, 
%      in which case the second column is used to specify the numerical label of 
%      the cluster.
%   IN.REF (string, default the AAL NIAK template) a file name of a 3D 
%      partition used as a reference to label the clusters. Note that the cluster
%      numerical values do not need to be 1, ... , N. The background is filled with 0s.
%   IN.LABELS_REF (string, necessary if REF is not the default, otherwise 
%      uses the AAL labels) same as LABELS_PART, but for REF.
%
% OUT (string) a .csv with a sumary table. Each row corresponds to a cluster:
%      column 1: string label (in IN.LABELS_PART)
%      column 2: numerical label (in IN.PART)
%      column 3: the volume (in mm3)
%      column 4: the label of the reference cluster with highest overlap with 
%         the cluster.
%      column 5: the percentage of overlap of the reference cluster inside the 
%         target cluster.
%      column 6: the percentage of overlap of the target clusters inside the 
%         reference cluster.
%      Each block of three sbsequent columns code for the second-most reference 
%      cluster with highest overlap, third-most reference cluster, etc. 
%      The max number of clusters can be set in OPT
%
% OPT.NB_REF (integer, default 5) the max number of reference clusters being listed. 
% OPT.FLAG_VERBOSE (boolean, default 1) if the flag is 1, then the function prints 
%   some infos during the processing.
% OPT.FLAG_TEST (boolean, default 0) if FLAG_TEST equals 1, the brick does not do 
%   anything but update the default values in IN, OUT and OPT.
%           
% _________________________________________________________________________
% OUTPUTS:
%
% IN, OUT, OPT: same as inputs but updated with default values.
%              
% _________________________________________________________________________
% COMMENTS:
%
% The target and reference partitions need to be in the same space, and sampled
% on the same grid !
% 
% An example of a label file would look like:
% Precentral_L  , 2001
% Precentral_R  , 2002
% Frontal_Sup_L , 2101
% etc 
% _________________________________________________________________________
% Copyright (c) Pierre Bellec
% Centre de recherche de l'institut de gériatrie de Montréal, 
% Department of Computer Science and Operations Research
% University of Montreal, Québec, Canada, 2010-2014
% See licensing information in the code.
% Maintainer : pierre.bellec@criugm.qc.ca
% Keywords : NIAK, table, clusters

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


%% Initialization and syntax checks
flag_gb_niak_fast_gb = true; % Only load the most important global variables for fast initialization
niak_gb_vars 

%% Syntax
if ~exist('in','var')||~exist('out','var')
    error('niak:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

if nargin < 3
    opt = struct;
end

%% Input files
in = psom_struct_defaults( in , ...
      {'part' , 'labels_part' , 'ref' , 'labels_ref' }, ...
      {NaN    , ''            , ''    , ''           });

if ~isempty(in.ref)&&isempty(in.labels_ref)
    error('Please specify FILES_IN.LABELS_REF to use FILES_IN.REF');
end

if isempty(in.ref)
    in.ref = [gb_niak_path_niak 'template' filesep 'roi_aal_3mm.mnc.gz'];
    in.labels_ref = [gb_niak_path_niak 'template' filesep 'labels_aal.csv'];
end

%% Output file
if ~ischar(out)
    error('FILES_OUT should be a string')
end

%% Options
opt = psom_struct_defaults( opt , ...
      {'nb_ref' , 'flag_verbose' , 'flag_test' }, ...
      {5        , true           , false       });

%% If the test flag is true, stop here !
if opt.flag_test == 1
    return
end

%% The core of the brick starts here

%% Read the target partition
[hdr,part] = niak_read_vol(in.part);
labels_part = niak_read_csv_cell(in.labels_part);
ind_part = unique(part(:));
ind_part = ind_part(ind_part~=0);
nb_part = length(ind_part);
if size(labels_part,2)>1
    [val,s_idx] = ismember(ind_part,labels_part(:,2));
    if any(~val)
        error('Some numerical labels in IN.PART could not be found in the second column of IN.LABELS_PART')
    end
    labels_part = labels_part(s_idx,1);
end
    
%% Read the reference partition
[hdr2,ref] = niak_read_vol(in.ref);
if any(size(part)~=size(ref))
    error('It looks like the target and reference clusters are not in the same space')
end
ind_ref = unique(ref(:));
ind_ref = ind_ref(ind_ref~=0);
nb_ref = length(ind_ref);
if size(labels_ref,2)>1
    [val,s_idx] = ismember(ind_ref,labels_ref(:,2));
    if any(~val)
        error('Some numerical labels in IN.PART could not be found in the second column of IN.LABELS_PART')
    end
    labels_ref = labels_ref(s_idx,1);
end

%% now make the table
tab = cell(nb_part+1,3+3*opt.nb_ref);

