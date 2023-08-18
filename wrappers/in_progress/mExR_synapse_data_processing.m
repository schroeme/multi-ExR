%%  WRAPPER for processing multi-ExR synapses data and converting to data table format%% 

% Accompanies the manuscript by Kang*, Schroeder* et al., 2023
% Last modified by Margaret Schroeder on 6/3/23

% Note: you can run this section-by-section or all at once.

clear all

% Enter directory where your .mat files (output of
% mExR_synapses_analysis_wrapper.m) here and change filenames appropriately
parentfolder = 'E:/Margaret/mExR/2023.03_synapses/cropped_rois/';
load([parentfolder 'S1-2_multiplex_imdatas_20230609.mat'])

%% Input protein and extract fov names

%edit your channel names here - these will get put into the data table
%for first mouse
% proteins = {'Bassoon';'SynGAP';'NR1';
%     'Shank3';'GluA1';'CaMKIIa';
%     'Cav2.1';'NR2B';'RIM1';
%     'PSD95';'Gephyrin';
%     'Vglut';'Vgat';
%     'RIMBP';
%     'GluA3';'Stargazin';
%     'GluA2';
%     };
proteins = {'SynGAP';'NR1';
    'RIM1';'Gephyrin';
    'GluA4';'IRSp53';
    'GluA1';'NR2B';
    'Homer1';'CaMKIIa';
    'Shank3';
    'Bassoon';
    'Erbb4';'Stargazin';
    'Elfn1';
    'PSD95';'Cav2.1';
    'Vglut';'GluA2';
    'GluA3';
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

alldatas = [all_im_datas];

%CHANGE THE FILENAME HERE
writetable(alldatas,[parentfolder 'allsynapsedata_S1-2_20230609.csv'])

