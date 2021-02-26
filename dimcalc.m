% Function Name: dimcalc.m

% Author 1: David Johnson
% Author 2: Timothy Secomb, PhD

% Version History:
%   1.0 - Initial creation

% Description
%   (1) Purpose: To generate dimensions based on parent-daughter
%       relationship in the seg tree structure with the capability to 
%       simulate vasoconstriction. 
%   (2) Algorithms or Techniques: N/A

% Input
%   seg - The tree structure for a network based on a number of terminal
%       vessels (capillaries) with multiple vessels coming out of one
%       parent vessel at the same node. This is later processed by
%       debrancher_art_tree to create branches off parent vessel by
%       breaking vessel into smaller segments to create new nodes within 
%       the parent vessel. 
%   rbart - The branching ratio specified for aterial circulation 
%       (presently, same as veins, but could be adapted to be different)
%   rbven - The branching ratio specified for venous circulation 
%       (presently, same as arteries, but could be adapted to be different)
%   dimrand - Determine vessel dimensions using a poisson distribution
%       based on murray's law (==0) or horsefield data (==1, not used).
%   squeezeord1 - Specifies smallest range of diameter reduction
%   squeezepct1 - Specifies % of original diameter of smallest range
%   squeezeord2 - Specifies largest range of diameter reduction
%   squeezepct2 - Specifies % of original diameter of largest range
%   murray - Used to specify indicies for dimension construction
%   level - Used to signify which level of network (1, 2, or 3) is being
%       simulated by current iteration
%   lseg - Used to specify last highest order of lower level network (not 
%       used for level 1 network, must be 1 higher than previous network's 
%       highest order)
%   ldiam - Used to speficy last highest vessel diameter of lower level 
%       network (not currently used, useful for debugging)
%   llen - Used to speficy last highest vessel length of lower level
%       network (not currently used, useful for debugging)
%   absvalsqueeze - Specifies whether constant or variable diameter
%       reduction is used
%   ncap - Number of capillaries supplied by the network

% Output 
%   segwdims - Tree structure with dimensions
%   xlowhigh - the low and high x values used between levels for plotting
%       in the cumulative resistance graph. 
%       (e.g. the xhigh of level 1 should be right next to the xlow of
%       level 2)


function [segwdims, xlowhigh] = dimcalc(seg, rbart, rbven, dimrand, squeezerange1, squeezepct1, squeezerange2, squeezepct2, murray, level, lseg, ldiam, llen, absvalsqueeze, ncap)  %maxord

% Assign length/diam of each segment
% REQUIRES ADDITIONAL FILE WITH HORSEFIELD DATA. Comment: not used due to
% questionable data and non-scaling with Murray's law; no longer necessary

% Assumes vein and artery rb is same, currently only arterial network is
%   assessed

%nseg = size(seg,1); 
%determine length of seg

xlowhigh = [-1,-1];

%Not currently used
%hors = readtable('C:\Users\david\Desktop\Secomb_Lab\Strahler_MATLAB\Strahler_MATLAB\Horsfield_Strahler_Data.csv');
% table format: hors = [order, hdiam, hleng]

%hors = table2array(hors);       %changes into array so find function works

if murray == 1      %Murray approximation, not currently used

    
% \\ 5/15
%     
%     caps = find(hors(:,1) == 0);
%     startD = hors(caps,2);
%     startL = hors(caps,3);
%     
%     hors = zeros(201,4);
%     for ik = 1:201
%         hors(ik,1) = 101-ik;
%     end
%     
%     hors(101,2) = startD;
%     hors(101,3) = startL;
%     
%     %caps = find(hors(:,1) == 0); %capillary row
%     ind = 0;
%     
% %     hors(caps,5) = hors(caps,2);
% %     hors(caps,6) = hors(caps,3);
%     
%     %diameter
%     for ord = 1:100
%         ind = ind+1;
%         hors(101 - ind,2) = power(rbart,abs(ord)/3)*hors(101,2);
%     end
%     ind = 0;
%     
%     for ord = -1:-1:-100
%         ind = ind+1;
%         hors(101 + ind,2) = power(rbven,abs(ord)/3)*hors(101,2);
%     end
%         
%     %length
%     ind = 0;
%     for ord = 1:100
%         ind = ind+1;
%         hors(101 - ind,3) = power(rbart,abs(ord)/3)*hors(101,3);
%     end
%     
%     ind = 0;
%     
%     for ord = -1:-1:-100
%         ind = ind+1;
%         hors(101 + ind,3) = power(rbven,abs(ord)/3)*hors(101,3);
%     end
%     

    m1_ind = 2;
    m2_ind = 3;
    
elseif murray == 0  %Murray scaling
    m1_ind = 2;
    m2_ind = 3;
else
    error('murray specification not 1 or 0')   
end


    for iseg = 1:size(seg(:,1))

        if (seg(iseg,1) > 0)
            rb = rbart;

            alphaD = 1/3;            % -- Murray's law on diameter
            alphaL = 1;              % -- NOT TRUE--> Murray's law doesnt apply to length
            
        elseif (seg(iseg,1) < 0)
            rb = rbven;
            
            alphaD = 1/3;            % -- Murray's law on diameter
            alphaL = 1;              % -- NOT TRUE --> Murray's law doesnt apply to length

        elseif (seg(iseg,1) == 0)
            rb = 1;
            % is this a bad assumption for branching ratio of capillaries?
            % also consider that we have two sides to capillaries (art and ven)
            % this takes away heterogeneity of capillaries

            alphaD = 0;            %no poisson for capillaries
            alphaL = 0; 
            
        end

        order = seg(iseg, 1);                  %Create Order var
        
        
        if dimrand == 1     %on
            diam = hors(find(hors(:,1) == seg(iseg,1)),m1_ind)*((seg(iseg,4)/rb)^alphaD);
            len = hors(find(hors(:,1) == seg(iseg,1)),m2_ind)*((seg(iseg,4)/rb)^alphaL);
            %Use poisson distr of number of branches diam = (hors estimate *
                % (Nb/rb)^alpha(for tuning) (nb = number of branches)
                
            seg(iseg,5) = diam;
            seg(iseg,6) = len;

        elseif dimrand == 0

            if level == 1
            
            diam = power(rbart,abs(seg(iseg,1))/3)*0.0006; % seg(iseg,1) = order
            len = power(rbart,abs(seg(iseg,1))/3)*0.0100;
% %                 Length between bifurcations should scale with diameter
% %                 Length for a given strahler order --> baseline length should increase when there is a higher branching ratio
% %                 Could we extract length:diameter ratio
            
            elseif level == 2
                
            % We used ordrel to reset the diameter at the beginning of each
            % level such that for level 1 -> level 2 the vessel is set to
            % 3808^1/3*(length or diameter) and scaled exactly based on the murray's law scaling
            % for the difference of order from lseg (the max order from L1) within a given level
            % We are doing this because just doing this by order is not
            % accurate enough and leads to systematic mistakes (using order
            % in the power function didnt work just with power alone)
            % If we scale by number of segments, then diameter or length
            % should scale by 3808^(1/3) for each level

            orderl = abs(seg(iseg,1))-lseg;    %This is the order within this level
                
            diam = power(rbart,orderl/3)*0.0006*(3808^(1/3)); % seg(iseg,1) = order
            len = power(rbart,orderl/3)*0.0100*(3808^(1/3));
            
% %             diam = power(rbart,seg(iseg,1)/3)*0.0006;
% %             len = power(rbart,seg(iseg,1)/3)*0.0100;
% %              
                
            elseif level == 3
              
            orderl = abs(seg(iseg,1))-lseg;    %This is the order within this level
                
            diam = power(rbart,orderl/3)*0.0006*(3808^(2/3)); % seg(iseg,1) = order
            len = power(rbart,orderl/3)*0.0100*(3808^(2/3));
           
% %             diam = power(rbart,seg(iseg,1)/3)*0.0006;
% %             len = power(rbart,seg(iseg,1)/3)*0.0100;
            
            end
            

%             Not currently used
%             diam = hors(find(hors(:,1) == seg(iseg,1)),m1_ind);
%             len = hors(find(hors(:,1) == seg(iseg,1)),m2_ind);

            if isempty(diam) == 1
               error('error') 
            end
            
            seg(iseg,5) = diam;
            seg(iseg,6) = len;
        
        else 
            error('dimrand setting unspecified')
        end
        
    end
    
    %Simulates conststriction -- squeezepct1 or squeezepct2 indicate % of
    %original diameter such that 100% is the original diameter and 90% is a
    %10% reduction
    if squeezepct1 ~= 100 || squeezepct2 ~=100
        
        if squeezerange1 == 1
            a = 0.0000;
        elseif squeezerange1 == 2
            a = 0.0020;
        elseif squeezerange1 == 3
            a = 0.0060;
        elseif squeezerange1 == 4
            a = 0.0180;
        elseif squeezerange1 == 5
            a = 0.0540;
        elseif squeezerange1 == 6
            a = 0.1620;
        elseif squeezerange1 == 7
            a = 0.4860;
        elseif squeezerange1 == 8
            a = 1.0000;
        end

        if squeezerange2 == 1
            b = 0.0020;
        elseif squeezerange2 == 2
            b = 0.0060;
        elseif squeezerange2 == 3
            b = 0.0180;
        elseif squeezerange2 == 4
            b = 0.0540;
        elseif squeezerange2 == 5
            b = 0.1620;
        elseif squeezerange2 == 6
            b = 0.4860;
        elseif squeezerange2 == 7
            b = 1.0000;
        elseif squeezerange2 == 8
            b = 10.0000;
        end
        
        if absvalsqueeze == 2
            squeezeslope = NaN;
        else
            %y variable
            squeezediffy = squeezepct2-squeezepct1;

            if squeezediffy ~=0
                error('squeeze percent mismatch')
            end

            %x variable is range of squeeze diameters
    %         maxdiam = max(seg(find(seg(:,5) < b),5));
    %         mindiam = min(seg(find(seg(:,5) > a),5)); 
            squeezediffx = b-a;  %maxdiam-mindiam;
            squeezeslope = squeezediffy/squeezediffx;
        end
        
        % Not used -- Vasoconstriction with sharp rise and fall without plateau 
        % this squeeze has a midpoint which takes on the 2nd % squeeze,
        % where the a and b values take on the 1st % squeeze
        if absvalsqueeze == 1
            n = 0;
            for iseg = 1:size(seg(:,1))
               if seg(iseg,5) >= a && seg(iseg,5) <= b
                   
                   %abs value graph where %2 is midpoint and %1 is at a and b
                   dd = abs(squeezepct1-squeezepct2)/(abs(b-a)/2);
                   ee = squeezediffx/2 + a;
                   ff = squeezepct2;
                   
                   squeezepctval = dd*abs(seg(iseg,5)-ee)+ff;
                   
                   %Debugging point
                   n = n+1;
                   squeezetrack(n)= squeezepctval;
                   %
                   
                   seg(iseg,5)= seg(iseg,5)*squeezepctval/100;
               end
            end
        end
        
        % Expected scenario: Variable vasoconstriction considering conducted response
        % best guess up ramp, plateau, down ramp
        if absvalsqueeze == 2 && squeezepct1 ~=100
           
            for iseg = 1:size(seg(:,1))
            
%                 if seg(iseg,5) >= 0.0050 && seg(iseg,5) < 0.0100
%                     squeezepctval = 100 - (30/log10(2))*(log10(seg(iseg,5))-log10(0.0050));
%                 elseif seg(iseg,5) >= 0.0100 && seg(iseg,5) < 0.0400
%                     squeezepctval = 70;
%                 elseif seg(iseg,5) >= 0.0400 && seg(iseg,5) < 0.0800
%                     squeezepctval = 70 + (30/log10(2))*(log10(seg(iseg,5))-log10(0.0400));
%                 else
%                     squeezepctval = 100;
%                 end

                if seg(iseg,5) >= 0.0033 && seg(iseg,5) < 0.0100
                    squeezepctval = 100 - ((100-squeezepct1)/log10(3))*(log10(seg(iseg,5))-log10(0.0033));
                elseif seg(iseg,5) >= 0.0100 && seg(iseg,5) < 0.0300
                    squeezepctval = squeezepct1;
                elseif seg(iseg,5) >= 0.0300 && seg(iseg,5) < 0.0900
                    squeezepctval = squeezepct1 + ((100-squeezepct1)/log10(3))*(log10(seg(iseg,5))-log10(0.0300));
                else
                    squeezepctval = 100;
                end

                seg(iseg,5)= seg(iseg,5)*squeezepctval/100;
           
            end
        end
        
        xlowfigure = 0;
        xhighfigure = 0;
        
        indexyyy = 0;
        indexzzz = 0;
        
        xlowseg = zeros(1,2);
        xhighseg = zeros(1,2);
        
%       max(seg(find(seg(:,5) < b)
        if absvalsqueeze == 0 || absvalsqueeze == 3 || absvalsqueeze == 5 
            for iseg = 1:size(seg(:,1))
               if seg(iseg,5) >= a && seg(iseg,5) < b

                   squeezepctval = squeezepct1 + squeezeslope*(seg(iseg,5)-a);  %mindiam);

                   seg(iseg,5)= seg(iseg,5)*squeezepctval/100;
               end
               
               if ncap == 3808
                   if seg(iseg,5) < a
                       indexyyy=indexyyy+1;
                       xlowseg(indexyyy) = seg(iseg,2);
                       %xlowfigure = xlowfigure+1;
                   end
                   if seg (iseg,5) < b
                       indexzzz=indexzzz+1;
                       xhighseg(indexzzz) = seg(iseg,2);
                       %xhighfigure = xhighfigure+1;
                   end
               end
               
            end
        end
        
        if ncap == 3808
            
            xlowfigure = size (unique(xlowseg),2);
            xhighfigure =  size (unique(xhighseg),2);
            xlowhigh = [xlowfigure,xhighfigure]*(3808^(level-1));
        else
            xlowhigh = [-1,-1];
        end
        
    end
    
segwdims = seg;

end