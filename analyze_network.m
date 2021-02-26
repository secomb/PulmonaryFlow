% Function Name: analyze_network.m

% Author 1: David Johnson
% Author 2: Timothy Secomb, PhD

% Version History:
%   1.0 - Initial creation

% Description
%   (1) Purpose: To create a tree structure based on nodal relations in
%       order to iteratively solve for nodal pressures using
%       solve_network.m
%   (2) Algorithms or Techniques: N/A

% Input
%   seg - The tree structure for a network based on a number of terminal
%       vessels (capillaries) with multiple vessels coming out of one
%       parent vessel that *has been modified* by debrancher_art_tree to 
%       create branches off parent vessel by breaking vessel into smaller 
%       segments to create new nodes within the parent vessel. 

% Output 
%   nod - The tree structure for a network based on nodal relations

function [ nod ]  = analyze_network (seg)

    maxnodes = max(max(seg(:,3:4)));
    nod = zeros(maxnodes,8);

    for inod = 1:maxnodes

        %why can't inod = nod name
        %inod = inod + 1;
        nod(inod,1) = inod;                 %inod is index

        %Nodseg
        nodinx = 0;

        segs = find(seg(:,3)==inod | seg(:,4)==inod);

        if isempty(segs)
                print('nodnod error, inx should be at least 1')
        else

            nodtype = length(segs);
            nod(inod,8) = nodtype;          %number of matches speficies type

            nodinx = 0;

            for iseg = segs'

                nodinx = nodinx + 1;

                %Assign connected segment as nodseg
                nod(inod,nodinx+1) = iseg;    

                %Assign connected nodes as nodnod
                if seg(iseg,3) ~= inod
                    %if 3rd collumn value is not same as iseg, then this is a connected node
                    nod(inod,4+nodinx) = seg(iseg,3);    %equiv to ista
                elseif seg(iseg,4) ~= inod
                    %if 4th collumn value is not same as iseg, then this is a connected node
                    nod(inod,4+nodinx) = seg(iseg,4);    %eqiv to iend
                end

            end

        end

    end

end