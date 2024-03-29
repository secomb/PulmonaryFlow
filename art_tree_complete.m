% Function Name: art_tree_complete.m

% Authors: David W. Johnson, Tuhin K. Roy and Timothy W. Secomb

% Version History:
%   1.0 - Initial creation

% Description
%   (1) Purpose: To create a strahler-based pulmonary arterial sub-network
%   structure in order to simulate the effects of vasoconstriction on
%   pulmonary vascular resistance. 
%   (2) Algorithms or Techniques: N/A

% Input
%   rb - Branching ratio of the network
%   ncap - Number of capillaries supplied by the network
%   dimrand - Determine vessel dimensions using a poisson distribution
%       based on murray's law (==0) or horsefield data (==1, not used).
%   squeezeord1 - Specifies smallest range of diameter reduction
%   squeezepct1 - Specifies % of original diameter of smallest range
%   squeezeord2 - Specifies largest range of diameter reduction
%   squeezepct2 - Specifies % of original diameter of largest range
%   murray - Used to specify indicies for dimension construction
%   seedon - Used to specify random number generator iteration (rng.m)
%   level - Used to signify which level of network (1, 2, or 3) is being
%       simulated by current iteration
%   lseg - Used to specify last highest order of lower level network (not 
%       used for level 1 network, must be 1 higher than previous network's 
%       highest order)
%   ldiam - Used to speficy last highest vessel diameter of lower level 
%       network (not currently used, useful for debugging)
%   llen - Used to speficy last highest vessel length of lower level
%       network (not currently used, useful for debugging)
%   ldratio_rb2 - The length:diameter ratio for all vessels in a network
%       created with branching ratio ==2 (deterministic)
%   absvalsqueeze - Specifies whether constant or variable diameter
%       reduction is used
%   randstatus - Used to specify random number generator method (rng.m)

% Output 
%   d_arttree_seg - Network structure generated by strahler.m, dimcalc.m,
%       and debrancher_art_tree.m
%   treestats - Used to communicate internal statistics (not currently
%       used)
%   network - A node-based network structure that indicates nodal
%       connection between segments, generated by analyze_network.m
%   netresis - The overall network resistance determined by solve_network.m
%   ldratio_rb2 - The length:diameter ratio for all vessels in a network
%       created with branching ratio ==2 (deterministic) -- Note: this is 
%       only an output generated when the branching ration ==2
%   xlowhigh - the low and high x values used between levels for plotting
%       in the cumulative resistance graph. 
%       (e.g. the xhigh of level 1 should be right next to the xlow of
%       level 2)

function [d_arttree_seg, treestats, network, netresis, ldratio_rb2, xlowhigh] = art_tree_complete (rb, ncap, dimrand, squeezeord1, squeezepct1, squeezeord2, squeezepct2, murray, seedon, level, lseg, ldiam, llen, ldratio_rb2, absvalsqueeze,randstatus)

%Obtain tree
[seg] = strahler(rb,ncap,seedon, level, lseg,randstatus);

if level == 1
    %Fix orders to include 
    for F = 1:length(seg(:,1))
        if (seg(F,1) == 1)                %Art-side Capillary order
            seg(F,1) = 0;
        elseif (seg(F,1) == -1)           %Ven-side Capillary order, not used in current project
            seg(F,1) = 0;
        elseif (seg(F,1) > 1)             %Arterial order
            seg(F,1) = seg(F,1) - 1;
        elseif (seg(F,1) < 1)             %Venous order, not used in current project
            seg(F,1) = seg(F,1) + 1;

        end
    end    
end
    
%Assign Parameters
[seg,xlowhigh] = dimcalc(seg, rb, rb, dimrand, squeezeord1, squeezepct1, squeezeord2, squeezepct2, murray, level, lseg, ldiam, llen, absvalsqueeze, ncap);

%Debranch art tree
d_arttree_seg = debrancher_art_tree(seg,level,lseg);

%%create nodnod and nodseg network for flow calculation
network = analyze_network(d_arttree_seg);

%%solve network
[~,~,~,netresis] = solve_network(network,d_arttree_seg,ncap,level,lseg);

% Statistical analysis for one newtork, debugging point
    % Check branching ratio to see if d^4/length * number of vessels = same accross orders
        %Uses not debranched seg (as opposed to d_arttree_seg)
    if ncap == 3808
        a = 0;
        for order = 0:max(seg(:,1))
               a = a+1;
               
               if level == 1
                   n = 3808^2;
               elseif level == 2
                   n = 3808;
               elseif level == 3
                   n = 1;
               end
               
               inx = find(seg(:,1) == order);
               avgdiam = mean(seg(inx,5));
               avglen = mean(seg(inx,6));
               check_ratio(a,1) = avgdiam^4/avglen*length(inx)*n;
               check_ratio(a,2) = order;
        end
    %     
    %     figure
    %     plot(check_ratio(:,2),check_ratio(:,1)) 
    %     
        if level > 1
            a = 1; %Debugging point
            
        end
        treestats = 1; %Not currently used
    end
        
end




