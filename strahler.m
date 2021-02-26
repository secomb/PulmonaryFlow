% Function Name: strahler.m

% Author 1: David Johnson
% Author 2: Timothy Secomb, PhD

% Version History:
%   1.0 - Initial creation

% Description
%   (1) Purpose: Create a strahler-like network where the number of
%       parent:number of daughter vessels at a single node is determined by
%       using a poisson distribution based on an average branching ratio.
%   (2) Algorithms or Techniques: N/A

% Input
%   rb - Branching ratio of the network
%   ncap - Number of capillaries supplied by the network
%   seedon - Used to specify random number generator iteration (rng.m)
%   level - Used to signify which level of network (1, 2, or 3) is being
%       simulated by current iteration
%   hi_prev_order - Used to specify last highest order of lower level 
%       network (not used for level 1 network, must be 1 higher than 
%       previous network's highest order) AKA lseg
%   randstatus - Used to specify random number generator method (rng.m)

% Output 
%   seg - The tree structure for a network based on a number of terminal
%       vessels (capillaries) with multiple vessels coming out of one
%       parent vessel at the same node. This is later processed by
%       debrancher_art_tree to create branches off parent vessel by
%       breaking vessel into smaller segments to create new nodes within 
%       the parent vessel. 

function [seg] = strahler(rb, ncap, seedon, level, hi_prev_order,randstatus)
    
    % Set initial number of segments
    nseg = zeros(100);
    nseg(1) = ncap;                 
    
    if seedon == 1                
    	%rng('default') -- same as rng(0)
        rng(randstatus)
    end
    
    % Create capillary segments
    iseg = 0;
    for j = 1:nseg(1)
        iseg = iseg + 1;                            % iseg accounts capillaries
        if level == 1 
            seg(iseg,1) = 1;                        % Capillaries are order 1 here (later they are order 0)
        elseif level == 2
            seg(iseg,1) = hi_prev_order;
        elseif level == 3
            seg(iseg,1) = hi_prev_order;
        else
            error('order not transferring to correct level')
        end
        seg(iseg,2) = j;                        % fromnode (CAPILLARY TO ARTERY/VEIN DIRECTION)
        seg(iseg,3) = 0;                        % tonode (TBA)
        seg(iseg,4) = 1;              %"TBD";   % No current method to express capillary branching
        
    end
    
    %index based on level
    if level == 1 
            i = 1;                        % Capillaries are order 1
        elseif level == 2
            i = hi_prev_order;
        elseif level == 3
            i = hi_prev_order;
    end
    
    while (find(seg(:,3)==0))         
        
        freeseg = find(seg(:,1) == i & seg(:,3) == 0); %assign freeseg to all seg(:,1) values where seg(:,3)=0 to free seg
        nleft = size(freeseg,1);
        k = 1;                                  % assigns first value of free seg later
        
        while (nleft > 0)                       % Cycle through vessels of order i
           
           %Determine br
           if (nleft == 3)                       %Ensure loop doesn't elongate network
               poissrb = nleft;
           elseif (nleft == 2)                   %Ensure loop doesn't elongate network 
               poissrb = nleft;
           else
               
               % New random branching ratio           
               mindex = poissrnd(rb - 2);        % Change loop to index random poisson with lambda = rb
               poissrb = mindex + 2;             % create poissrb variable to record each branching ratio
                    
           end 
           
           if (nleft-poissrb == 1)               %Precents having rb where 1 vessel is left
                   poissrb = nleft;         
                  
           elseif (poissrb > nleft)              %Precent rb greater than nleft
                   poissrb = nleft;
                   
           end
           
           if (nleft > 1)
                % Construct new segment of order i+1 AND from node number iseg+1
                nseg(i+1) = nseg(i+1) + 1;          % increase number of vessels in new order
                iseg = iseg + 1;                    % iseg has already accounted capillaries
    
                seg(iseg,1) = i + 1;                % Order of NEW segment
                seg(iseg,2) = iseg;                 % fromnode = NEW segment node
                seg(iseg,3) = 0;                    % tonode (TBA)
                seg(iseg,4) = poissrb;
    
                %Assign current iseg to tonode for segment 1 out of # branches (rb) 
                seg(freeseg(k),3) = iseg;           % tonode = NEW node
                nleft = nleft - 1;                  % one vessel was used
                k = k + 1;                          % next free seg
    
                %Assign current iseg to tonode for segment 2 out of # branches (rb) 
                seg(freeseg(k),3) = iseg;           % tonode = NEW node
                nleft = nleft - 1;                  % one vessel was used
                k = k + 1;                          % next free seg
           
                % Add additional segments where rb > 1
                if (nleft >= 1)
    
                    for rbleft = 1:(poissrb - 2)
    
                        %Assign to node for segment 2
                        seg(freeseg(k),3) = iseg;           % tonode = NEW node
                        nleft = nleft - 1;                  % one vessel was used
                        k = k + 1;                          % next free seg
    
                    end 
                end      
           end
            
            
            if (nleft == 1)
                
                
                
                %deal with situation where one segment is left
                if (seg(iseg,1) ~= max(seg(:,1)))         % Indicates not final order
                    seg(freeseg(k),3) = iseg;                   % tonode = NEW node
                    nleft = nleft - 1;                          % one vessel was used
                    seg(iseg,4) = poissrb + 1;                     % FIXES branching ratio to accomodate +1 branch
                        
                    
                    if (nleft ~= 0)
                        print('nleft error')
                        
                    end
                
                %final branch needs to be used and while loop ends ends    
                elseif (seg(iseg,1) == max(seg(:,1)))     % Indicates final order
                    seg(iseg,3) = iseg+1;                       % this assigns final branch a final to-node
                    break                                       % this breaks while loop
                    
                end            
            end
    
            % Create new set of free segs when nleft = 0
            if (nleft == 0)
               
                i = i + 1;                                      % increases order by 1           
                freeseg = find(seg(:,1) == i & seg(:,3) == 0);  % create new set of free segs for next order
                nleft = size(freeseg,1);
                k = 1;                                          % assigns first value of free seg later    
                
            end
            
        end
    
    end
    
    %Debugging point
    %g = digraph(seg(:,2),seg(:,3))
    %plot(digraph(seg(:,2),seg(:,3)))            % ncap = 50 tends to plot nicely but higher is worse
    
    
    
    nnode = size(seg,1) + 1;                        % total number of nodes [output]
    imax = i - 1;                                   % records imax minues last 
    
    artseg = seg;
    totnseg = size(seg,1);
    
    for i = 1:totnseg                            
        if (seg(i,2) <= nnode)
            seg(i,2) = nnode + 1 - seg(i,2);
        end
        if (seg(i,3) <= nnode)
            seg(i,3) = nnode + 1 - seg(i,3);
        end
    end
    % Flip arterial side of tree to put input segment on top
    % this is literally just putting artseg in reverse order

    %Flips to and from node
    for i = 1:totnseg
        artseg(i,1) = seg(totnseg+1-i,1);
        artseg(i,2) = seg(totnseg+1-i,3);
        artseg(i,3) = seg(totnseg+1-i,2);
        artseg(i,4) = seg(totnseg+1-i,4);
    end
    
    for i = 1:totnseg
        seg(i,:) = artseg(i,:);
    end
end

