% Function Name: solve_network.m

% Authors: David W. Johnson, Tuhin K. Roy and Timothy W. Secomb

% Version History:
%   1.0 - Initial creation

% Description
%   (1) Purpose: To determine a solution of network pressures of internal
%       nodes based on arbitrary boundary node pressures and specified
%       dimensions within parent-daugter vessel relationship. The network
%       pressures are used to determine network resistance and may be 
%       compared to observed changes in pulmonary resistance due to 
%       diameter changes (e.g. Hypoxic pulmonary vasoconstriction or 
%       vascular remodeling)
%   (2) Algorithms or Techniques: Overrelaxation method

% Input
%   nod - The tree structure for a network based on nodal relations
%   seg - The tree structure for a network based on a number of terminal
%       vessels (capillaries) with multiple vessels coming out of one
%       parent vessel that *has been modified* by debrancher_art_tree to 
%       create branches off parent vessel by breaking vessel into smaller 
%       segments to create new nodes within the parent vessel. 
%   ncap - Number of capillaries supplied by the network
%   level - Used to signify which level of network (1, 2, or 3) is being
%       simulated by current iteration
%   lseg - Used to specify last highest order of lower level network (not 
%       used for level 1 network, must be 1 higher than previous network's 
%       highest order)

% Output 
%   nodpress - The node pressure matrix that is iteratively adjusted,
%       references nod matrix
%   flowsum - The sum of flows for a given node, this is iteratively
%       ajusted to be 0 by through overrelaxation adjustment of node 
%       pressures
%   mp_pres - The midpoint pressure between two nodes, assumed at the
%       center of a vessel segment
%   netresis - The net resistance of a whole aterial network (including 
%       the contribution of lower level netwrok structure resistance)

function [nodpress, flowsum,mp_pres,netresis] = solve_network(nod,seg,ncap,level,lseg)

omega = 1;               %if this gets too high may not converve
nnod = length(nod(:,1)); %max number of nodes
nodsegm = max(nod(:,8)); %max segments out of nodes
wk = zeros(nodsegm,nnod); % sets up coefficient matrix from 1:nodsegm, 1:nnod %working coefficient matrix --> in C++ its dmatrix(1,nodsegm,1,nnod)
%wk starts off as conductance matrix and become coefficient matrix

maxerr = 1e6; %set large initial maxerr
niter = 0;  %number of iterations
nitmax = 10^4;
tolerance = 10^-5;

% Conductance Calculation % from Poiseuille R = 8 nu L / pi r^4
% Assumes constant viscosity (and hct)
% facfp = pi * 1333 (dyn/cm^2 per mmHg) * 128 (p's law) * 0.001 (P viscosity nu) * 60 (sec/min) * 10^6 (um^3/nL)
% (dyn)(cm^-2)(poise = dyn cm/s^2)(60s/1min)(um^3/nL)
facfp = pi()/128; %*1333/128/0.001*60/1e6; 
visc = 0.03; %3 cP = 0.03 P
nseg = length(seg(:,1));
cond = zeros(nseg,1);

%Boundary Hydrostatic pressures - defines type 1 nodes can have inlet or outlet conditions
%(West's Respiratory Physiology, 10th edition, p42 figure 4.1) --> 
% mean pulmonary arterial wedge pressure should be ~9 
%(this pressure is negative in lung)
P_in = 15*1333;     %dyne/cm^2  %Mean pulmonary arterial pressure is 15mmHg but 25 during systole, 8 during diastole       
P_out = 5*1333;    %ranges from 8-12mmHg from arterial to venous end

errnode = 0;

%output array
nodpress = zeros(length(nod(:,1)),1);


press1 = 0;

for iseg = 1:length(seg(:,1))
    cond(iseg) = facfp*(seg(iseg,5)^4)/seg(iseg,6)/visc;
    
    % Set conductance to be super high (resistance contribution is low) for 
    % level 2 or 3 lseg vessels to reduce the effects of double counting 
    % vessels and smooth jumps between levels
    if level > 1
        if seg(iseg,8) == lseg
            cond(iseg) = cond(iseg)*10000;
        end
    end
end


%Set up coefficients
for inod = 1:nnod

    % conditional nod type statements
    if nod(inod,8) == 1
        wk(1,inod) = 0;
    
    elseif nod(inod,8) > 1   
        condsum = 0;                            %Conductance sum for a given node
        for i = 1:nod(inod,8)
            iseg = abs(nod(inod,1+i));
            condsum = condsum + cond(iseg);     %Where is this cond coming from?
            wk(i,inod) = cond(iseg);
        end
        
        for i = 1:nod(inod,8)
           wk(i,inod) = wk(i,inod)/condsum;
        end
    end
end

bcprfl = zeros(ncap+1,2);
%setup inodbc matrix
for inod = 1:nnod
    if nod(inod,8) == 1
        
        bcprfl(inod,1) = inod;
        
        onlyseg = nod(inod,2);
        
        if level == 1
            i = 0;
        elseif level >1
            i = lseg;
        end
        
        if seg(onlyseg,8) == i
            bcprfl(inod,2) = P_out;
            wk(1,inod) = P_out;
            if nnod == 2 && inod == 1   %This is for a network of 1 segment and 2 nodes ONLY
                bcprfl(inod,2) = P_in;
                wk(1,inod) = P_in;
            end
            
        else
            bcprfl(inod,2) = P_in;
            wk(1,inod) = P_in;
        end
        
    end
end


%Set Pressure/Flow Boundary Nodes -> temporarily set nodtype of pressure nodes to -1
for inodbc = 1:length(bcprfl(:,1))
   
    inod = bcprfl(inodbc,1);                         
    if nod(inodbc,8) == 1
        
        nodpress(inod) = bcprfl(inodbc,2);  
        nod(inod,8) = -1;
        
        
    end
    
     %%%else <<< do we even need this if else statement if bctype is always
     %%%boundary?
        %%%wk(1,inod) = bcprfl(inodbc)/cond(nod(inod,1)); %% should P_in just be 1 instead?
    
        %still confused about bcpfl
        %I interpreted it to mean that P_in = 1 and P_out = 0 and 
        %all other working matrix values are a proportion of Pin
        %Why is wk only referenced in 1st row for else statement?
        
        %bcpressureflow is an array that gives boundary pressures
        % need boundary conditions array with pout and pin 
        % just need 2vector array with bcnod values and bc pressure values [1,0] arbitrary

end



%Iterative process for pressures
while maxerr > tolerance && niter <= nitmax

    maxerr = 0;                     %maxerr is updated every loop where > 0
    for inod = 1:nnod
        
        if nod(inod,8) == -1
        
            press1 = wk(1,inod) - nodpress(nod(inod,5));
            %press1 = omega*(nodpress(nod(inod,5)) + wk(1,inod));  % omega = overrelaxation method
        
        elseif nod(inod,8) >= 2
            
            pcondsum = 0; %sum of wk(1)*p1 + wk(2)*p2 + wk(3)*p3
            for i = 1:abs(nod(inod,8))
                %there's something wrong below bc nodpress is calling
                %attached nodes to inod
                %this is actually fine
                pcondsum = pcondsum + wk(i,inod)*nodpress(nod(inod,4 + i)); %% --> this is flow
                
                %%%*** why is nodpress --> how is nodpres defined for
                %%%non-boundary condition nodes
                %%% is nodpress created as loop cycles? and it doesn't
                %%% change for boundary conditions?
                
                press1 = omega*(pcondsum - nodpress(inod));    
            end
        %end
        
        %if nod(inod,8) >= 1
        
            %Updating nod pressure
            nodpress(inod) = nodpress(inod) + press1; %<< This is output of function
            
            %evaluate/store max err
            if abs(press1) > maxerr
                maxerr = abs(press1);
                errnode = inod;
            end
        end
    end
    
    niter = niter + 1;
    if niter == nitmax
           
        msg1 = strcat('excessive number of loops (', num2str(niter), ') max error of (', num2str(maxerr), ')', 'Check Omega Value');
        error(msg1)
        break
        
    end
    
end

for inod = 1:nnod 
    if nod(inod,8) == -1
        nod(inod,8) = 1;
    end
end
     
%Flow calc > should see ~zero in all non-boundary nodes
flowsum = zeros(inod,1);
for inod = 1:nnod
        
        %Obtain all vessels leading to node
        seg1 = nod(inod,2);
        seg2 = nod(inod,3);
        seg3 = nod(inod,4);
        nod1 = nod(inod,5);
        nod2 = nod(inod,6);
        nod3 = nod(inod,7);
        
        if nod(inod,8) == 1
            
            f1 = cond(seg1)*(nodpress(nod1)-nodpress(inod));
            flowsum(inod) = f1;
            
        elseif nod(inod,8) == 3
      
            %Assume all flows are into node
            f1 = cond(seg1)*(nodpress(nod1)-nodpress(inod));
            f2 = cond(seg2)*(nodpress(nod2)-nodpress(inod));
            f3 = cond(seg3)*(nodpress(nod3)-nodpress(inod));
            flowsum(inod) = f1 + f2 + f3;        
        
        end
    
end

%calculate midpoint pressures and flow in each segment
for iseg = 1:nseg

    fromnod = seg(iseg,3);
    tonod = seg(iseg,4);
    mp_pres(iseg) = (nodpress(fromnod)+nodpress(tonod))/2;
    flow(iseg) = cond(iseg).*(nodpress(fromnod,1)-nodpress(tonod,1));        %?????
    
end

%Effective resistance is deltaP/flow in inlet segment (flow is conserved)
% Think of resistors in parallel: 1/R1+1/R1 = 1/Reff simplified to Reff = R1/2
localnetworkresis = (P_in-P_out)/flow(1);
totlungcaps = 5.5248*10^8/(0.01);     %From secomb and roy (2014) dividing total length by the assumed capillary length of 100um 
totnumnetworks = totlungcaps/ncap; %total number of networks of this size ncap in lung

% If statements are wrong and unncessary
if level == 2
    totnumnetworks = totnumnetworks/(3808);
elseif level == 3
    totnumnetworks = totnumnetworks/(3808^2);
end
    
netresis = localnetworkresis/totnumnetworks; %this is the net resistance from order 0 to #; dyne/sec/cm^5
        %%% totnumnetworks should be totlungcaps divided by (3808^2 = ncap) for second hop 
        %%% should be 3808^3 for next
        %%% resistance at end of first hop --> should be first hop + resistance caulculated at second hop

netresis;

end