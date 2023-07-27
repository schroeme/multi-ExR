function res = run_nanocolumn_enrichment_refAlign(params,fov)
%update later with documentation

exp_factor = params.exp_factor;   %expansion factor, adjust as needed
parentdir= params.parentdir;
targetdir = params.targetdir;
nChannels = params.nChannels;
rmax = params.rmax;
prepost=params.prepost; %assign pre or post-synaptic identity to channel. 0 = pre, 1 = post
pixel_size=params.pixel_size;
step = params.step;
channels_rounds = params.channels_rounds;

res = [];       %array to keep the results
ch01_syn_files = dir([parentdir fov '*round001*ch01*.tif']);
nsynapses = length(ch01_syn_files); %number of synapses for this ROI

for sss = 1:nsynapses %loop through all synapses
    syn_name_splits = split(ch01_syn_files(sss).name,'_');
%     syn_name = [syn_name_splits{5} '_' syn_name_splits{6}];
    syn_name = syn_name_splits{4};
    syn_files = dir([parentdir fov '*' syn_name]);
    
    for ii = 1:nChannels%loop through all channels
        for jj = 1:nChannels%do this in a pairwise manner
            
            ch1str=channels_rounds{ii};
            ch1splits = split(ch1str,'-');
            ch1roundstr = ['0' ch1splits{1}];
            ch1str = ch1splits{2};
           
            ch1files = dir([parentdir fov '*round*' ch1roundstr '_ch0' ch1str '*.tif']);
            file1 = fullfile(ch1files(sss).name);

            ch2str=channels_rounds{jj};
            ch2splits = split(ch2str,'-');
            ch2roundstr = ['0' ch2splits{1}];
            ch2str = ch2splits{2};
           
            ch2files = dir([parentdir fov '*round*' ch2roundstr '_ch0' ch2str '*.tif']);
            file2 = fullfile(ch2files(sss).name);

            test1 = loadtiff([parentdir file1]);
            test1 = double(test1);

            test2 = loadtiff([parentdir file2]);
            test2 = double(test2);
            
            %Subdivide pixel size (physical pixel size of camera chip divided by
            %magnification of objective) and step size of stack in order to achieve
            %isometric sub-voxels. 
            im1 = expand(test1,[2,2,3]); %2,2,3 are divisors that give equal pixel values for x,y,z
            im2 = expand(test2,[2,2,3]); 
            pixel = pixel_size/exp_factor; %Pixel size defined here, in nm

            if (prepost(ii) == prepost(jj))%(strcmp(ch1roundstr,ch2roundstr)) && (prepost(ii) == prepost(jj)) %if we're in the same round, and same pre-post identity, no xyz shift
                xyzshift = [0,0,0];
            else 
                xyzshift = [];
            end
        
            %Shift clusters to be compared so that they are overlapped. If clusters
            %are expected to be colocalized and occupy the same space (and thus require
            %no shift), set xyzshift = [0,0,0]
            distance = params.distance;
            flag = params.flag;
            
            [Raa,Rab,Rba,Rbb] = get_enrichment_3dMatrix_final(im1, im2, pixel, rmax, xyzshift, distance, step, flag);
            restemp = [1,sss,ii,jj,Rab'];
%             abtemp = [sss,ii,jj,Rab'];
%             batemp = [sss,ii,jj,Rba'];
%             bbtemp = [sss,ii,jj,Rbb'];
            
            res = [res;restemp];
%             resaa = [resaa;aatemp];    %A density along distance to peak of A
%             resab = [resab;abtemp];    %A density along distance to peak of B
%             resba = [resba;batemp];    %B density along distance to peak of A
%             resbb = [resbb;bbtemp];    %B density along distance to peak of B
        end
    end

end
save([targetdir fov '_AtoB_Enrichment_refAlign.txt'], 'res', '-ascii', '-tabs');
% save('AtoBpeak.txt', 'resab', '-ascii', '-tabs');
% save('BtoApeak.txt', 'resba', '-ascii', '-tabs');
% save('BtoBpeak.txt', 'resbb', '-ascii', '-tabs');

%Column 1: ROI id
%Column 2: Synapse number in that ROI
%Column 3: Channel 1
%Column 4: Channel 2
%Column 5-end: Enrichment results
end