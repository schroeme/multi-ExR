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

res = [];       %array to keep the results
ch01_syn_files = dir([parentdir fov '*round001*ch01*.tif']);
nsynapses = length(ch01_syn_files); %number of synapses for this ROI

for sss = 1:nsynapses %loop through all synapses
    syn_name_splits = split(ch01_syn_files(sss).name,'_');
    syn_name = [syn_name_splits{5} '_' syn_name_splits{6}];
    syn_files = dir([parentdir fov '*' syn_name]);
    
    for ii = 1:nChannels%loop through all channels
        for jj = 1:nChannels%do this in a pairwise manner
            file1 = fullfile(syn_files(ii).name);
            splits1 = split(file1,'_');
            round1 = splits1{2};

            file2 = fullfile(syn_files(jj).name);
            splits2 = split(file2,'_');
            round2 = splits2{2};

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

            if (strcmp(round1,round2)) && (prepost(ii) == prepost(jj)) %if we're in the same round, and same pre-post identity, no xyz shift
                xyzshift = [0,0,0];
            elseif (strcmp(round1,round2)) && (prepost(ii) ~= prepost(jj)) %if we're in the same round, and diff pre-post identity, yes xyshift
                xyzshift = [];
            elseif (strcmp(round1,round2)==0) && (prepost(ii) == prepost(jj)) %if we're in different rounds, and same pre/post identity, shift only to reflect reg error
%                 ref1files = dir([parentdir fov '*' round1 '*ch04*' syn_name_splits{5} '_' syn_name_splits{6}]);
%                 ref1file = ref1files(1).name;
%                 ref2files = dir([parentdir fov '*' round2 '*ch04*' syn_name_splits{5} '_' syn_name_splits{6}]);
%                 ref2file = ref2files(1).name;
% 
%                 ref1 = loadtiff([parentdir ref1file]);
%                 ref1 = double(ref1);
%                 ref2 = loadtiff([parentdir ref2file]);
%                 ref2 = double(ref2);
%                 im1 = expand(test1,[2,2,3]); %2,2,3 are divisors that give equal pixel values for x,y,z
%                 im2 = expand(test2,[2,2,3]); 
%                 pixel = pixel_size/exp_factor; %Pixel size defined here, in nm
%                 
%                 xyzshift = [];
%                 distance = [10,800]; %minimum and maximum shift radius, based on registration error
%                 flag = 0;
%                 xyzshift = calculate_xyzshift(ref1, ref2, pixel, rmax, xyzshift, distance,flag);
                xyzshift = [];
            elseif  (strcmp(round1,round2)==0) && (prepost(ii) ~= prepost(jj)) %allow any xyzshift
                %later we can figure out if a double-shift is needed
                xyzshift = [];
            end
        
            %Shift clusters to be compared so that they are overlapped. If clusters
            %are expected to be colocalized and occupy the same space (and thus require
            %no shift), set xyzshift = [0,0,0]
            distance = [10,800]; %Min and max distance allowed for xyz shift, in nm
            flag = 0;
            
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