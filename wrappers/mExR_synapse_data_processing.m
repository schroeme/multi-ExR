%%  WRAPPER for processing multi-ExR synapses data and converting to data table format%% 

% Accompanies the manuscript by Kang*, Schroeder* et al., 2023
% Last modified by Margaret Schroeder on 12/14/2022

% Note: you can run this section-by-section or all at once.

clear all

% Enter directory where your .mat files (output of
% mExR_synapses_analysis_wrapper.m) here and change filenames appropriately
parentfolder = 'E:/Margaret/mExR/2022.08_synapses/synapses_cropped/';
load([parentfolder 'multiplex_imdatas_20221119.mat'])
load([parentfolder 'multiplex_syndatas_20221119.mat'])
load([parentfolder 'multiplex_pwdatas_20221119.mat'])

%% Input protein and extract fov names

%edit your channel names here - these will get put into the data table
proteins = {'Bassoon';'SynGAP';'NR1';'Homerref1';
    'Shank3';'GluA1';'CaMKII';'Homerref2';
    'Cav2.1';'NR2B';'RIM1';'Homerref3';
    'PSD95';'Gephyrin';'Homerref4';
    'Vglut';'Vgat';'Homerref5';
    'mGluR5';'PSD95R6';'Homerref6';
    'RIMBP';'Adam22';'Homerref7';
    'GluA3';'Stargazin';'Homerref8';
    'Synapsin1';'GABAB';'Homerref9';
    'Synaptophysin';'NR1R10';'Homerref10';
    'GluA4';'GluA2';'Homerref11'
    };
nproteins = length(proteins);

for p1 = 1:length(proteins)
    for p2 = p1+1:length(proteins)
        pair = [proteins{p2} '-' proteins{p1}];
        pairs{p1,p2} = pair;
    end
end
pairs =  pairs(~cellfun('isempty',pairs));

for im = 1:length(imdatas)
    name = imdatas(im).filename;
    splits = split(name,"\");
    fovs{im} = splits{1};
end

%fovs = fovs([1 2]); %manually exclude fovs that were not well aligned
%imdatas = imdatas([1 2]);
%syndatas = syndatas([1 2]);
%pwdatas = pwdatas([1 2]);

%% collapse pairwise matrix data into rows

fields = fieldnames(pwdatas(1).pwdata);
A = reshape([1:nproteins^2],nproteins,nproteins);
tri_inds = triu(A,1);
tri_inds = tri_inds(tri_inds>0);
tri_inds = sort(tri_inds);

clear labeled_data namedpairs
for im = 1:length(fovs) %loop images
    for ss = 1:size(pwdatas(im).pwdata,2) %loop through synapses
        labeled_datas = [];
        for ff = 1:length(fields) %loop through fields
            datamat = getfield(pwdatas(im).pwdata(ss),fields{ff});
            tridata = datamat(tri_inds);
            varNames = {'variable',[fovs{im} '_syn' num2str(ss)]};
            for pp = 1:length(pairs) %relabel the pair name
                namedpairs{pp,1} = [fields{ff} '_' pairs{pp}];
            end
            labeled_data = table(namedpairs,tridata,'VariableNames',varNames);
            labeled_datas = [labeled_datas; labeled_data];
        end
        if ss == 1
            combined_datas = labeled_datas;
        else
            combined_datas = join(combined_datas,labeled_datas);
        end
    end
    if im == 1
        all_pw_datas = combined_datas;
    else
        all_pw_datas = join(all_pw_datas,combined_datas);
    end
end

%% collapse syndata matrices data into rows

fields = fieldnames(syndatas(1).syndata);
A = reshape([1:nproteins^2],nproteins,nproteins);
tri_inds = triu(A,1);
tri_inds = tri_inds(tri_inds>0);
tri_inds = sort(tri_inds);

clear labeled_data namedpairs
for im = 1:length(fovs) %loop images
    for ss = 1:size(syndatas(im).syndata,2) %loop through synapses
        labeled_datas = [];
        for ff = 1:length(fields) %loop through fields
            datamat = getfield(syndatas(im).syndata(ss),fields{ff});
            tridata = datamat(tri_inds);
            varNames = {'variable',[fovs{im} '_syn' num2str(ss)]};
            for pp = 1:length(pairs) %relabel the pair name
                namedpairs{pp,1} = [fields{ff} '_' pairs{pp}];
            end
            labeled_data = table(namedpairs,tridata,'VariableNames',varNames);
            labeled_datas = [labeled_datas; labeled_data];
        end
        if ss == 1
            combined_datas = labeled_datas;
        else
            combined_datas = join(combined_datas,labeled_datas);
        end
    end
    if im == 1
        all_syn_datas = combined_datas;
    else
        all_syn_datas = join(all_syn_datas,combined_datas);
    end
end

%% reformat imdata into a table

fields = fieldnames(imdatas(1).imdata);

clear labeled_data namedpairs data
for im = 1:length(fovs) %loop images
    for ss = 1:size(imdatas(im).imdata.mean_puncta_vol,1) %loop through synapses
        labeled_datas = [];
        for ff = 1:length(fields) %loop through fields
            data = getfield(imdatas(im).imdata,fields{ff});
            varNames = {'variable',[fovs{im} '_syn' num2str(ss)]};
            namedprots = {};
            if contains(fields{ff},'size')
                data = data(ss,1);
                namedprots{1,1} = fields{ff};
            else
                data = data(ss,:)';
                for pp = 1:length(proteins) %relabel the pair name
                    namedprots{pp,1} = [fields{ff} '_' proteins{pp}];
                end
            end
            labeled_data = table(namedprots,data,'VariableNames',varNames);
            labeled_datas = [labeled_datas; labeled_data];
        end
        if ss == 1
            combined_datas = labeled_datas;
        else
            combined_datas = join(combined_datas,labeled_datas);
        end
    end
    if im == 1
        all_im_datas = combined_datas;
    else
        all_im_datas = join(all_im_datas,combined_datas);
    end
end

%% combine all data tables and write to csv

alldatas = [all_im_datas; all_syn_datas; all_pw_datas];

%CHANGE THE FILENAME HERE
writetable(alldatas,[parentfolder 'allsynapsedata_ROI4_20221120.csv'])

