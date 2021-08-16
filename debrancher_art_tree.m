% Function Name: debrancher_art_tree.m

% Authors: David W. Johnson, Tuhin K. Roy and Timothy W. Secomb

% Version History:
%   1.0 - Initial creation

% Description
%   (1) Purpose: To create intermediate branch points to covert seg tree
%       structure to represent purely bifurcating tree with parent vessels 
%       supplying 2 or more daughter vessesl
%   (2) Algorithms or Techniques: N/A

% Input
%   seg - The tree structure for a network based on a number of terminal
%       vessels (capillaries) with multiple vessels coming out of one
%       parent vessel at the same node. This is later processed by
%       debrancher_art_tree to create branches off parent vessel by
%       breaking vessel into smaller segments to create new nodes within 
%       the parent vessel. 
%   level - Used to signify which level of network (1, 2, or 3) is being
%       simulated by current iteration
%   lseg - Used to specify last highest order of lower level network (not 
%       used for level 1 network, must be 1 higher than previous network's 
%       highest order)

% Output 
%   arttree - tree structure that has "debranched" vessels by creating
%       intermediate branch points in parent vessels such that each node 
%       has no greater than 3 vessels connected

function [arttree] = debrancher_art_tree(seg,level,lseg)
% New Vessels Array Size(includes split same-order vessels)

if level >1
    for k = 1:length(seg(:,1))
        if seg(k,1) > 0 
            seg(k,1) = seg(k,1) - lseg;
            seg(k,7) = 1;                   %indicates arterial 
%         elseif seg(k,1) < 0 
%             seg(k,1) = seg(k,1) + lseg;
%             seg(k,7) = -1;                  %indicates venous, not used
        else
            error('venous tree should not be in seg')
        end
    end
end


dseg_number = 0;
for i = 1:length(seg(:,1))
    if (seg(i,4) <= 2)                                  %Branching ratio less than 2 contributes no additional vessels
        dseg_number = dseg_number + 1;
    elseif (seg(i,4) > 2)                               %Contributes additional new vessels when rb>2
        dseg_number = dseg_number + (seg(i,4)-1);
    end
end
dseg = zeros(dseg_number,9);                            %Create new array of size of new debranched network

renumber = 10^ceil(log10(dseg_number));

% Indicies
nseg = 0;                                                %index for segment number on dseg
nnode = 1;                                               %index for to-node number (starts at 1) on dseg
rbindex = 0;                                             %index for rb if >2
%oldseg is an index for seg not dseg
ncap = 0;


% For Loop renumbering and debranching vessels 
for oldseg = 1:length(seg(:,1))
    
    rb = seg(oldseg,4);
    
    
% Loop enters debrancher if rb > 2 for retrograde vessel     
    if rb >2    
        
        rbindex = rb;
        branchcounter = 0;
        
        while rbindex >= 2
                    
            nseg = nseg + 1;
            dseg(nseg,1) = nseg;                        %Assign New seg #
            dseg(nseg,2) = 0;                           %type = 0 as default
            
            if branchcounter == 0                       %from-node assigment
                dseg(nseg,3) = seg(oldseg,2);           %From-node --> Same but different number scheme
                                                        %First node will be = 1 which doesn't matter 
            elseif branchcounter > 0
                dseg(nseg,3) = dseg(nseg-1,4);          %from node of seg is referenced
                
            end
            
            nnode = nnode + 1;                                        
            dseg(nseg,4) = nnode + renumber;    
                                                        %Assign to-node ***rb dependence
                                                        %To-node is of one less order than from-node
                                                        %ASSUMPTION: There is no order skipping among vessels
            
            %Accounting in seg to rename from-nodes
            for inx = find(seg(:,2)==seg(oldseg,3))
                if rbindex > 2                          
                    
                    if inx == 0
                        break
                    end
                    
                    a = isempty(inx);
                    
                    if a == 1
                        break
                    end
                    
                    seg(inx(1),2) = dseg(nseg,4);       %Assign anterograde vessel from old # scheme to apropriate new from node
                    
                    rbindex = rbindex - 1;
                    
                end
            end
            
            dseg(nseg,5) = seg(oldseg,5);               %Transfer diam
                    
            dseg(nseg,6) = seg(oldseg,6)./(rb-1);       %*** rb dependence if rb>2
                                                        %Transfer length 

            dseg(nseg,7) = seg(oldseg,4);               %Transfer branching ratio
            dseg(nseg,8) = seg(oldseg,1);               %Transfer order  
            
            if level >1
                dseg(nseg,9) = seg(oldseg,7);            %Arterial/venous indication
            end

                
            %Creating NEW anterograde segments of SAME ORDER/DIMENSIONS -- divide old seg
            %ONLY if rb index indicates >2 anterograde lower order branches remain
            while rbindex > 1
                nseg = nseg + 1;
                nnode = nnode + 1;
                
                dseg(nseg,:) = dseg(nseg-1,:);                      %Create new anterograde seg of same
                dseg(nseg,1) = nseg;                                %Update seg #
                dseg(nseg,3) = dseg(nseg-1,4);                      %to-node of retrograde = from-node of anterograde
                dseg(nseg,4) = nnode + renumber;                       
                
                for inx = find(seg(:,2)==seg(oldseg,3))
                                      
                    if rbindex > 2
                    
                        seg(inx(1),2) = nnode + renumber;       %Assign anterograde vessel from old # scheme to apropriate new from node
                    
                        rbindex = rbindex - 1;
                        branchcounter = branchcounter + 1;
                        
                    elseif rbindex == 2
                    
                        dseg(nseg,4) = nnode + renumber;
                        seg(inx,2) = nnode + renumber;
                                                        %Assgn both remaining anterograde lower order vessels
                        
                        rbindex = rbindex - 2;
                        branchcounter = branchcounter + 1;
                        
                    end
                end
            end
            
            while rbindex == 1
                'error rbindex should not equal 1'
                
            end
          
        end
            
% Loop enters normal numbering scheme if rb = 2 for retrograde vessel
    elseif rb <= 2
            
        nseg = nseg + 1;
        dseg(nseg,1) = nseg;                        %Assign New seg #
        dseg(nseg,2) = 0;                         %type = 0 as default

        dseg(nseg,3) = seg(oldseg,2); 
                                                    %From-node --> Same but different number scheme
                                                    %The from-node should be changed to new numbering scheme however, first node will be = 1 which doesn't matter 

        nnode = nnode + 1;                                        
        dseg(nseg,4) = nnode + renumber;  %***rb dependence
                                                %To-node is of one less order than from-node
                                                %ASSUMPTION: There is no order skipping among vessels

                for inx = find(seg(:,2)==seg(oldseg,3))

                    seg(inx,2) = dseg(nseg,4);

                end

        dseg(nseg,5) = seg(oldseg,5);             %Transfer diam

        dseg(nseg,6) = seg(oldseg,6);         %*** rb dependence
                                                %Transfer length 

        dseg(nseg,7) = seg(oldseg,4);             %Transfer branching ratio
        dseg(nseg,8) = seg(oldseg,1);             %Transfer order
        if level >1
                dseg(nseg,9) = seg(oldseg,7);            %Arterial/venous indication
        end

    end

end

dseg(1,3) = dseg(1,3) + renumber;
dseg(:,3:4) = dseg(:,3:4) - renumber;

arttree = dseg;

if level >1
   arttree(:,8) = arttree(:,8) + lseg; 
end

end