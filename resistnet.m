% Function Name: resistnet.m

% Authors: David W. Johnson, Tuhin K. Roy and Timothy W. Secomb

% Version History:
%   1.0 - Initial creation

% Description
%   (1) Purpose: Used to develop figures based on use of art_tree_complete
%       to establish 3 separate levels representing the pulmonary arterial
%       tree to simulate the effects of vasoconstriction (diameter 
%       reduction) on pulmonary resistance
%   (2) Algorithms or Techniques: N/A

% Input
%   rblow - Lowest branching ratio to use in simulations 
%   rbhigh - Highest branching ratio to use in simulations
%   highestpower - Specified the largest network of capillaries simualated
%   estHPV - Used to estimate HPV if ==100 (Not used in gen_net.m) 
%   squeezeord1 - Specifies smallest range of diameter reduction
%   squeezepct1 - Specifies % of original diameter of smallest range
%   squeezeord2 - Specifies largest range of diameter reduction
%   squeezepct2 - Specifies % of original diameter of largest range
%   n - Specifies plotting setting (n == 1 in gen_net.m)
%   seedon - Used to specify random number generator iteration (rng.m)
%   deltR - Used the track change in resistance between levels
%   absvalsqueeze - Specifies whether constant or variable diameter
%       reduction is used
%   ldratio_correction - Used to correct lenth:diameter ratio based on
%       bifurcating network to scale withe Murray's Law
%   randstatus - Used to specify random number generator method (rng.m)

% Output 
%   deltares - Tracks change in resistance between levels 
%   L3_resistance - Indicates resistance of whole network 
%   ldratio_correction - Used to correct length:diameter ratio to scalp
%       resistance appropriately with Murray's Law

function [deltares, L3_resistance, ldratio_correction] = resistnet(rblow, rbhigh, highestpower, estHPV, squeezeord1, squeezepct1, squeezeord2, squeezepct2,n,seedon, deltR, absvalsqueeze,ldratio_correction,randstatus)
%rblow and rbhigh can be non-integer branching ratios
%highest power should be total number of log10(3808) where 3808 is first of 3 levels in creating overall network

boxdims = zeros(3,2,8)-0.5;

%create outer rb, diameter
if n == 1
    if absvalsqueeze == 3 || absvalsqueeze == 2 || absvalsqueeze == 5
    else
        figure('name','rblines','DefaultAxesFontSize',15)
    end
end

%must specify nn=0
if n > 0
    subplot(1,1,n)
end 

rbind=0;
totlungcaps = 5.5248*10^8/(0.01);


if absvalsqueeze == 3 || absvalsqueeze == 2
    totallines = 100;
elseif absvalsqueeze == 5
    totallines = 100; 
else
    totallines = (rbhigh-rblow)/0.5+1;
end
colmap = hsv(totallines);

ldratio_rb2 = 1;

for i = rblow:0.5:rbhigh

    rb = i;
    rbind = rbind+1;
    lseg = 0;

    for level = 1:3
    
        if level == 1
            jind = 0;
            ldiam = 0.0006;
            llen = 0.0100;
        end
        
        if level == 2
            lseg = maxOrd1;
        elseif level == 3
            lseg = maxOrd2;
        end
        
        for j = 0:(highestpower)/50:highestpower 
            % Note: you can change the sampling rate of the graphs by
            % changing the denominator
            ncap = round(10^j);
            jind = jind+1;
            [art_tree_seg,~,~,resistance, ldratio_rb2, xlowhigh] = art_tree_complete (rb,ncap,0,squeezeord1, squeezepct1, squeezeord2, squeezepct2,1,seedon, level, lseg, ldiam, llen, ldratio_rb2, absvalsqueeze,randstatus); %Resistance, random off, murray on
            
            ldiam = max(art_tree_seg(:,5));
            llen = max(art_tree_seg(:,6));
            
            if level == 1
                rbmatrix(jind,1,rbind) = resistance;
                maxR1 = resistance;
                maxOrd1 = art_tree_seg(1,8);
            elseif level == 2
                rbmatrix(jind,1,rbind) = resistance + maxR1;
                maxR2 = resistance + maxR1;
                maxOrd2 = art_tree_seg(1,8);
            elseif level == 3
                rbmatrix(jind,1,rbind) = resistance + maxR2;
            end
                       
            rbmatrix(jind,2,rbind) = rb;
            
            if level == 1
                rbmatrix(jind,3,rbind) = log10(ncap);
            elseif level == 2
                rbmatrix(jind,3,rbind) = log10(ncap)+log10(3808);
            elseif level == 3
                rbmatrix(jind,3,rbind) = log10(ncap)+log10(3808)*2;
            end
            
            rbmatrix(jind,4,rbind) = art_tree_seg(1,5); %diameter
        end
        
        
                    % Length/diameter correction factor
                    % this factor is used to scale resistance in order to include the
                    % expected scaling of length with increased branching ratio
                    if level == 1
                        if squeezepct1 == 100 && squeezepct2 == 100   
                           inx = find(art_tree_seg(:,8)~=0);              
                           ldratio_rbany = mean(art_tree_seg(inx,6)./art_tree_seg(inx,5)); 

                           if rb == 2
                               ldratio_rb2 = ldratio_rbany;
                               ldratio_correction(rbind) = 1;
                            elseif rb ~= 2          
                               ldratio_correction(rbind) = ldratio_rb2/ldratio_rbany;
                           end
                        end
                    end

                    if j == highestpower
                        boxdims(level,:,rbind) = xlowhigh;
                    end
    end
    
    
    
                    % correction to allow length/diameter ratio to scale 
                    % with rb length/diameter ratio should scale with rb 
                    % between bifurcations based on Murray's emprical 
                    % finding, not trifurcations so this corrects the 
                    % impact of debranching which must allow a constant 
                    % l/d ratio --> this essentially increases resistance 
                    % as would occur (poiseueille) if length were to 
                    % increase as a result of increased branching ratio
                    % see dimcalc for scaling factor                
                    rbmatrix(:,1,rbind) = rbmatrix(:,1,rbind)*ldratio_correction(rbind);
    
    if level == 3
        
        L3_resistance(rbind,1) = rb;
        L3_resistance(rbind,2) = squeezeord1;
        L3_resistance(rbind,3) = squeezeord2;
        L3_resistance(rbind,4) = squeezepct1;
        L3_resistance(rbind,5) = squeezepct2;
        L3_resistance(rbind,6) = rbmatrix(jind,1,rbind); 
        
    end
      
    %creates new figure if plotting multiple graphs
    if n > 0
    else
        if i == 1 
            figure
        end
    end
    
    hold on
    
    if absvalsqueeze == 3 || absvalsqueeze == 1 || absvalsqueeze == 2
        if rbind == 3
        hold on
        rblines(rbind,1) = plot(rbmatrix(:,3,rbind),rbmatrix(:,1,rbind),'Color',colmap(squeezepct1,:));
        rblabels(rbind,1) = i;
        hold off
        end
    elseif absvalsqueeze == 5 
        if rbind == 3
            hold on
            rblines(rbind,1) = plot(rbmatrix(:,3,rbind),rbmatrix(:,1,rbind),'Color',colmap(squeezepct1,:));
            rblabels(rbind,1) = i;
            hold off
            ylim([0 200])
            xlim([0,12.25])
            
            
            %Creating rectangle to represent constriction range
            if boxdims(3,1,rbind) == 1
                
               if boxdims(2,1,rbind) == 1
                   xlow = log10(boxdims(1,1,rbind));
                   
               else 
                   xlow = log10(boxdims(2,1,rbind)) + log10(3808);
               end
            elseif boxdims(3,1,rbind) == -1
                xlow = -1;   
            else
                xlow = log10(boxdims(3,1,rbind)) + log10(3808)*2;
            end
            
            if boxdims(3,2,rbind) == 1
                
               if boxdims(2,2,rbind) == 1
                   xhigh = log10(boxdims(1,2,rbind));
                   
               else 
                   xhigh = log10(boxdims(2,2,rbind)) + log10(3808);
               end
            elseif boxdims(3,2,rbind) == -1
                xhigh = -1;
            else
                xhigh = log10(boxdims(3,2,rbind)) + log10(3808)*2;
            end
            
            if xlow ~= -1
            hold on    
            end
            
%             %Mild resistance elevation from hypoxia (less effects of increased CO) in humans (baseline multiplied by factor)
%             %Create HPV Reference from 100% diam from Talbot
%             %https://journals.physiology.org/doi/pdf/10.1152/japplphysiol.00903.2004
%             if squeezepct1 == 100
%                 hpvx = [10.4];
%                 %56 is est HPV elevation see below
%                 raw = [1.1 1.8];
%                 factor = raw./1.1;
%                 
%                 hpvresistance = L3_resistance(rbind,6)*factor;
%                 hpvy = [hpvresistance];
%                 hold on
%                 plot(hpvx,hpvy,'kp')
%                 hold off
%             
%             
%             %Spectrum of resistance elevation factor in canines from zero to complete hypoxia
%             %Create HPV Reference from 100% Diam from Marshall 1994
%             %https://www-sciencedirect-com.ezproxy1.library.arizona.edu/science/article/pii/0034568794901295
%                 hpvx = [10.4];
%                 %Order of raw is 100mmHg O2, 50, 0
%                 raw = [1.2 1.8 3.3];
%                 factor = raw./1.2;
%                 
%                 %56 is est HPV elevation see below
%                 hpvresistance = L3_resistance(rbind,6)*factor;
%                 hpvy = [hpvresistance];
%                 hold on
%                 plot(hpvx,hpvy,'kd')
%                 hold off
%             
% %             %Create Elevation of PVR with exercise (same as normal)
% %             %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3204788/
% %             %PVR remains constant at 1.2 Wood units (mmHg/L/min) w/ exercise
% %             %1.2 Wood units = 1.2 (80) Dynes sec cm^-5 = 96
% %                 hpvx = [10.4];
% %                 %56 is est HPV elevation see below
% %                 hpvresistance = 96;
% %                 hpvy = [hpvresistance];
% %                 hold on
% %                 plot(hpvx,hpvy,'gs')
% %                 hold off
% %             end
            
        end
    else
        rblines(rbind,1) = plot(rbmatrix(:,3,rbind),rbmatrix(:,1,rbind),'Color',colmap(rbind,:));
        rblabels(rbind,1) = i;
    end
    
    xlabel('Log_1_0 (number of terminal vessels)','FontSize',20)
    ylabel('Resistance (dynes-seconds-cm^-^5)','FontSize',20)
    set(0,'defaultAxesFontSize',15)
    
end

k = find(rblabels(:,1)==0);
    rblabels(k,1) = NaN;

    
if n == 0
    
    hold on
    legend(rblines(1:rbind,1),'Branching Ratio '+string(rblabels(:,1)),'fontsize',20)
    xlabel('log10 (Capillary Number)','fontsize',15); ylabel('Resistance (dynes-sec-cm^-5)','fontsize',20);
    title(strcat('Pulmonary Resistance v. Pulmonary Capillaries',' Sqeeze Range:',string(squeezeord1),' %Diam:',string(squeezepct1)),'fontsize',20);
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;

else

    if squeezeord1 == 1
        a = '0 um';
    elseif squeezeord1 == 2
        a = '25 um';
        a = '20 um';
    elseif squeezeord1 == 3
        a = '50 um';
        a = '60 um';
    elseif squeezeord1 == 4
        a = '100 um';
        a = '180 um';
    elseif squeezeord1 == 5
        a = '200 um';
        a = '540 um';
    elseif squeezeord1 == 6
        a = '400 um';
        a = '1620um';
    elseif squeezeord1 == 7
        a = '800 um';
        a = '4860 um';
    elseif squeezeord1 == 8
        a = '1 cm';
    end

    if squeezeord2 == 1
        b = '25 um';
        b = '20 um';
    elseif squeezeord2 == 2
        b = '50 um';
        b = '60 um';
    elseif squeezeord2 == 3
        b = '100 um';
        b = '180 um';
    elseif squeezeord2 == 4
        b = '200 um';
        b = '540 um';
    elseif squeezeord2 == 5
        b = '400 um';
        b = '1620um';
    elseif squeezeord2 == 6
        b = '800 mm';
        b = '4860 um';
    elseif squeezeord2 == 7
        b = '1 cm';
    elseif squeezeord2 == 8
        b = '10 cm';
    end
    
    if absvalsqueeze == 2
        title(strcat('Best Guess Rb:',num2str(rb),'-',num2str(squeezepct1)))
    else
        title(strcat('[',a,']:[',b,']','-[',string(squeezepct1),'%]:[',string(squeezepct2),'%]'),'fontsize',20)
    end
end

hold off


if estHPV == 100
    %Expected TPR box --> bound by Q = 4-6 L/s, Pin-Pout = 5mmHg P/Q=R
    %Pulmonary arterial wedge pressure varies from 2-15 but is typically 9
    NormalTPR = [log10(totlungcaps) 20; log10(totlungcaps) 130];

    %Dcided to just use values from literature 20-130
    %%[log(totlungcaps) 2*80/(4) ; log(totlungcaps) 2*80/(5);log(totlungcaps) 2*80/(8);...
    %            log(totlungcaps) 15*80/(4) ; log(totlungcaps) 15*80/(5);log(totlungcaps) 15*80/(8)];

    hold on
    plot(NormalTPR(:,1),NormalTPR(:,2),'DisplayName','Typical Range', 'LineWidth', 8, 'color', [0,0,0]); %1mmHg min/L == 80 dyn sec/cm^5
    hold off

    %Corrected Expected Change in TPR from Hypoxia (minus effects of increasing
    %CO to blunt increase of resistance take TPRc (c means corrected) from 
    %Table 6 B (right) resistance change with HPV corrected to remove effects of CO TPR = 0.7mmHg/L min (plot digitizer) 
    %From Talbot et Al (2004) https://journals.physiology.org/doi/pdf/10.1152/japplphysiol.00903.2004
    HPVTPRval = 0.7*80; %%1mmHg min/L == 80 dyn sec/cm^5
    HPVTPR(:,1) = NormalTPR(:,1)+.5;
    HPVTPR(:,2) = NormalTPR(:,2) + HPVTPRval;

    hold on
    plot(HPVTPR(:,1),HPVTPR(:,2),'DisplayName','HPV Range', 'LineWidth', 8, 'color', [.5,.5,.5]);
    hold off
end

if n == -1
    deltares =  2*rbmatrix(jind,1,:)-rbmatrix(1,1,:);
else 
    deltares = deltR;
end
 

if squeezepct1 == 100
 txt = ['\leftarrow' '  '  num2str(100-squeezepct1) '%'];  
else 
 txt = ['\leftarrow' num2str(100-squeezepct1) '%'];
end

if absvalsqueeze == 2
text(11,rbmatrix(jind,1,rbind)+1,txt,'FontSize',13, 'HorizontalAlignment','left')    
else
text(10.8,rbmatrix(jind,1,rbind)+1,txt,'FontSize',13, 'HorizontalAlignment','left')
end

end