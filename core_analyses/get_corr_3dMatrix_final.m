function out = get_corr_3dMatrix_final(A1, B1, pixel, rmax, xyzshift, distance,flag)


% function [xyzshift, d, 0, g] = get_crosscorr_3dMatrix(A, B, pixel, rmax, xyzshift,distance,flag)
% calculates crosscorrelation function between two 3d Matrix
%       For two clusters supposed to be overlaped, set xyzshift=[0,0,0];
%       For two clusters not overlaped, set xyzshift=[] if the best
%       overlaping vector is unknow. Use distance=[min, max] to limit the
%       |xyzshift|. The function will find the vector xyzshift to get the
%       best overlap between min(A,Amax/mv) and min(B,Bmax/mv). 
%
% INPUTS
% A = matrix of channel A 
% B = matrix of channel  B
% pixel = voxel size of the matrix
% rmax = maximum r value to correlate in units of pixels. default is 100;
% xyzshift = xyz vector to overlap the two clusters.
% distance = [min, max], distance variation allowed for xyzshift. for RIM
%     and PSD95, we use [40,160]. This also incorporate any potential error
%     for dual channel registration.
% flag = display flag.  insert 1 to display errorbar(r, g, dg) after
%   computation. 
%
% OUTPUTS
% xyzshift = [dx, dy, dz]
% d = |xyzshift|
% r = radius values
%
% NOTE: G(r=0) is just the dot product of the image.  For display purposes,
% G(r=0) is set to zero in the 3D autocorrelation output.  g(r=0) [g(1)]
% retains the proper value.
%
% Updated 01.26.10 by Sarah Veatch.
% Modified 02.18.14 by Aihui Tang.
% Last updated 10.21.17 by Aihui Tang.

if nargin~=7, flag = 0; end  % flag for display
n = 5;
k7 = zeros(n,n,n);
for i = 1:n
    for j = 1:n
        for k = 1:n
            if sqrt((i-(n+1)/2)^2+(j-(n+1)/2)^2+(k-(n+1)/2)^2)<=n/2
                k7(i,j,k) = 1;
            end
        end
    end
end
k7 = k7/sum(k7(:));

bka = median(A1(:));
bkb = median(B1(:));
k = [0,1,0;1,1,1;0,1,0]/5;
A1 = max(0,convn(A1-bka,k,'same')); 
B1 = max(0,convn(B1-bkb,k,'same'));

A2 = max(0,A1-max(A1(:))/3); B2 = max(0,B1-max(B1(:))/3); 
tha = mean(A2(find(A2>0))); thb = mean(B2(find(B2>0))); 
A2 = max(0,A1-tha); B2 = max(0,B1-thb); 
A2 = convn(A2,k7,'same'); B2 = convn(B2,k7,'same');
maskA = double(A2>0); maskB = double(B2>0); 
A2 = A1.*maskA; B2 = B1.*maskB;
A3 = min(A2,max(A2(:))/6); B3 = min(B2,max(B2(:))/6); 

Na = sum(sum(sum(A2.*maskA)));  % number of particles within channel 1
Nb = sum(sum(sum(B2.*maskB)));  % number of particles within channel 2
 Va = sum(sum(sum(maskA)));      % volume of maskA
 Vb = sum(sum(sum(maskB)));      % volume of maskB
 Va = sum(sum(sum(A3)));      % volume of maskA
 Vb = sum(sum(sum(B3)));      % volume of maskB
NP = real(fftshift(ifftn(fftn(maskA,size(A2)+rmax).*conj(fftn(maskB,size(B2)+rmax)))));
NP = real(fftshift(ifftn(fftn(A3,size(A2)+rmax).*conj(fftn(B3,size(B2)+rmax)))));
G1 = Va*Vb/Na/Nb*real(fftshift(ifftn(fftn(A2.*maskA,size(A2)+rmax).*conj(fftn(B2.*maskB,size(B2)+rmax)))))./NP; % 2D G(r) with proper normalization
L0=size(G1);L=floor(L0/2+1);
NPmask = max(0,max(max(max(NP)))/4-NP);
G1(find(NPmask))=0;             %set non-relavent points to 0

md = round(160/pixel);   %detect peak of correlation only in area with d<160 nm to center
if isempty(xyzshift) 
    G0 = real(fftshift(ifftn(fftn(A3.*maskA,size(A3)+rmax).*conj(fftn(B3.*maskB,size(B3)+rmax))))); 
    G0(find(NPmask))=0;
    G = G0*0; 
    pG = find(G==0);
    [x,y,z]=ind2sub(size(G), pG);
    d = pdist2([x, y, z],L);
    temp = find(d > distance(1)/pixel & d < distance(2)/pixel);
    G(pG(temp)) = 1;
    G = G.*G0; 
    Gm = max(max(max(G)));
    [mx, my, mz] = ind2sub(size(G), find(G==Gm));         %find the peak of correlation
    xyzshift = [mx(1)-L(1),my(1)-L(2),mz(1)-L(3)]; 
end
    %xyzshift=xyzshift.*(pixel/ps);
    mxyz = floor(xyzshift + L0/2+1);
    
g = nan(1,rmax+1);
if mxyz(1)-rmax>1 & mxyz(1)+rmax<size(G1,1) & mxyz(2)-rmax>1 & mxyz(2)+rmax<size(G1,2) & mxyz(3)-rmax>1 & mxyz(3)+rmax<size(G1,3)
    G = G1(mxyz(1)-rmax:mxyz(1)+rmax, mxyz(2)-rmax:mxyz(2)+rmax, mxyz(3)-rmax:mxyz(3)+rmax);

xvals = ones(rmax*2+1,rmax*2+1,rmax*2+1); yvals =xvals; zvals =xvals;
for i=1:2*rmax+1
    xvals(i,:,:)=i-rmax-1;
    yvals(:,i,:)=i-rmax-1;
    zvals(:,:,i)=i-rmax-1;%map to x positions with center x=0
end
r = sqrt(xvals.^2 + yvals.^2 + zvals.^2);

Ar = reshape(r,1, (2*rmax+1)^3);
Avals = reshape(G, 1, (2*rmax+1)^3);
temp = find(Avals ~=0);
Ar = Ar(temp); Avals = Avals(temp);
[rr,ind] = sort(Ar);                         % sort by r values
vv = Avals(ind);                             % reindex g
r = 0:floor(max(rr));                        % the radii you want to extract
[n bin] = histc(rr, r-.5);                   % bin by radius

dg = g;
for j = 1:rmax+1;                            % now get averages
    m = bin==j;
    n2 = sum(m);                             % the number of pixels in that bin
    if n2==0, vals(j)=0; er(j)=0;            % if no bins, no data
    else
        g(j) = sum(m.*vv)/n2;               % the average G values in this bin
        dg(j) = sqrt(sum(m.*(vv-g(j)).^2))/n2; % the variance of the mean
    end
end

r = 0:rmax;

%end

G(rmax+1, rmax+1, rmax+1) = 0;

if flag,
    r = 0:rmax;
    figure('Color', 'white');
    errorbar(r(1:length(r)), g(1:length(r)), dg(1:length(r)));
    axis tight
end

end
xyzshift;
out=[xyzshift Na/Va Nb/Vb g]; 